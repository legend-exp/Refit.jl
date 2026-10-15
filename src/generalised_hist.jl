"""
    GeneralisedHistogram(spectrum; grid...)
    GeneralisedHistogram(hists; grid...)

A histogram whose counts depend on parameters. It holds the histogram simulated at
every point of a grid in the parameters, and interpolates linearly between the grid
points, in each parameter and separately in each bin. [`get_weights`](@ref) and
[`get_histogram`](@ref) return the counts at given values of the parameters.

The keyword arguments give the name and the grid of each parameter, an evenly spaced
range of at least two values. All the histograms must have the same binning. They are
given in one of two ways:

- `spectrum` is a function that returns the histogram for the values of the parameters,
  given as a named tuple:
  ```julia
  gh = GeneralisedHistogram(p -> simulate(p.fccd, p.tl); fccd = 0.6:0.2:1.4, tl = 0.0:0.1:0.5)
  ```
- `hists` is an array of histograms with one dimension per parameter, in the order of
  the keyword arguments: `hists[i, j]` is the histogram at the `i`-th value of the first
  parameter and the `j`-th value of the second. Only the size of the array is checked,
  so two parameters with grids of the same length can be swapped by mistake. The
  function form above avoids this.

Each evaluation sums the histograms of the ``2^N`` corners of the grid cell around the
parameter values, so its time doubles with each parameter.

!!! warning
    The interpolation is an approximation. It suits a parameter that changes the counts
    in each bin smoothly. It does not suit a parameter that moves a peak, such as the
    energy scale: between two grid points the peak appears twice, at half the height.
    Check that the grid is fine enough for each use.

# Fields
* `edges`: the bin edges, in the units of the histograms (keV for a spectrum)
* `grid`: a named tuple with the grid values of each parameter
* `counts`: an array with the counts per bin at each grid point, `counts[i, j]`
  belongs to the `i`-th value of the first parameter and the `j`-th of the second
"""
struct GeneralisedHistogram{N,G<:NamedTuple}
    edges::Vector{Float64}
    grid::G
    counts::Array{Vector{Float64},N}
end

function GeneralisedHistogram(spectrum::Function; grid...)
    grid = NamedTuple(grid)
    # the histogram at every combination of grid values, in the order of the keywords
    hists = map(Iterators.product(Tuple(grid)...)) do values
        spectrum(NamedTuple{keys(grid)}(values))
    end
    return GeneralisedHistogram(hists; grid...)
end

function GeneralisedHistogram(
    hists::AbstractArray{<:Histogram{<:Real,1},N};
    grid...,
) where {N}
    grid = NamedTuple(grid)
    length(grid) >= 1 || throw(ArgumentError("a generalised histogram needs a parameter"))
    length(grid) == N || throw(
        ArgumentError("hists has $N dimensions but there are $(length(grid)) parameters"),
    )
    for (name, values) in pairs(grid)
        values isa AbstractRange ||
            throw(ArgumentError("the grid of $name must be a range, such as 0.5:0.1:2.0"))
        length(values) >= 2 && step(values) > 0 ||
            throw(ArgumentError("the grid of $name needs at least two increasing values"))
    end
    size(hists) == map(length, Tuple(grid)) || throw(
        DimensionMismatch(
            "hists has size $(size(hists)), the grid $(map(length, Tuple(grid)))",
        ),
    )

    edges = collect(Float64, first(hists).edges[1])
    all(h -> h.edges[1] == edges, hists) ||
        throw(ArgumentError("all histograms must have the same binning"))
    counts = map(h -> collect(Float64, h.weights), hists)
    return GeneralisedHistogram(edges, grid, counts)
end

# The cell of the grid `values` that contains `value`: the index of its lower grid
# value, and the position of `value` inside the cell, from 0 at the lower to 1 at the
# upper grid value
function _grid_cell(values::AbstractRange, value::Real, name::Symbol)
    first(values) <= value <= last(values) || throw(
        DomainError(
            value,
            "$name is outside the grid, from $(first(values)) to $(last(values))",
        ),
    )
    # the distance from the first grid value in grid steps. Rounding can put it a
    # hair outside the grid, so the cell index is kept inside
    x = (value - first(values)) / step(values)
    low = clamp(floor(Int, x), 0, length(values) - 2)
    return low + 1, x - low
end

"""
    get_weights(gh::GeneralisedHistogram, params::NamedTuple)

Counts per bin of `gh` at the parameter values in `params`, interpolated between the
grid points (see [`GeneralisedHistogram`](@ref)). `params` must contain every
parameter of the grid, in any order, and may contain others, which are ignored. A
value outside the grid is an error.
"""
function get_weights(gh::GeneralisedHistogram{N}, params::NamedTuple) where {N}
    names = keys(gh.grid)
    cells = ntuple(d -> _grid_cell(gh.grid[d], params[names[d]], names[d]), Val(N))
    low = map(first, cells)  # index of the lower grid value of the cell
    t = map(last, cells)     # position inside the cell, from 0 to 1

    # the element type of the parameters, so that derivatives (ForwardDiff) work
    weights = zeros(promote_type(Float64, map(typeof, t)...), length(gh.edges) - 1)
    # add the spectrum of each corner of the cell: 0 is the lower and 1 the upper grid
    # value along a parameter
    for corner in Iterators.product(ntuple(_ -> (0, 1), Val(N))...)
        w = prod(ntuple(d -> corner[d] == 1 ? t[d] : 1 - t[d], Val(N)))
        counts = gh.counts[(low .+ corner)...]
        for i in eachindex(weights)
            weights[i] += w * counts[i]
        end
    end
    return weights
end

"""
    get_histogram(gh::GeneralisedHistogram, params::NamedTuple)

The histogram of `gh` at the parameter values in `params`, see [`get_weights`](@ref).
"""
get_histogram(gh::GeneralisedHistogram, params::NamedTuple) =
    Histogram((gh.edges,), get_weights(gh, params))

"""
    rebin_hist(gh::GeneralisedHistogram, binning)

Rebin the spectrum at every grid point of `gh`, with an integer number of bins to
merge or a vector of new bin edges, as for a `Histogram`. Rebinning and interpolation
are both sums of counts, so rebinning before interpolating gives the same result as
interpolating before rebinning.
"""
function rebin_hist(gh::GeneralisedHistogram, binning::Union{Int,AbstractVector{<:Real}})
    hists = map(c -> rebin_hist(Histogram((gh.edges,), c), binning), gh.counts)
    return GeneralisedHistogram(hists; gh.grid...)
end

"""
    add_histograms(ghs::Vector{<:GeneralisedHistogram})

Sum the counts of generalised histograms with the same binning and the same grid,
grid point by grid point.
"""
function add_histograms(ghs::Vector{<:GeneralisedHistogram})
    all(gh -> gh.edges == ghs[1].edges, ghs) || error("Histogram binning mismatch")
    all(gh -> gh.grid == ghs[1].grid, ghs) || error("Grid mismatch")
    return GeneralisedHistogram(ghs[1].edges, ghs[1].grid, sum(gh.counts for gh in ghs))
end
