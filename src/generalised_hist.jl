"""
    GeneralisedHistogram(hists; grid...)

A simulated energy spectrum that depends on 1 to 5 parameters, for example the
thickness of a dead layer or the position of a source. It holds the spectrum
simulated at every point of a grid in these parameters and interpolates between
them.

`hists` is an array of 1D histograms with one dimension per parameter: `hists[i, j]`
is the spectrum simulated at the `i`-th value of the first parameter and the `j`-th
value of the second. All histograms must have the same binning. The keyword
arguments give the name and the grid values of each parameter, in the order of the
dimensions of `hists`. Each grid is an evenly spaced range of at least two values,
for example

```julia
gh = GeneralisedHistogram(hists; fccd = 0.5:0.1:2.0, alpha_frac = 0.0:0.2:1.0)
```

for a grid of 16 FCCD values (in mm) times 6 values of `alpha_frac`.

[`get_weights`](@ref) returns the interpolated counts per bin at any point inside
the grid. In every bin the counts are interpolated linearly in each parameter
(multilinear interpolation, the extension of linear interpolation to several
parameters). With ``N`` parameters the point lies in a cell of the grid with
``2^N`` corners. Along parameter ``d`` it lies a fraction ``t_d`` (from 0 to 1) of
the way from the lower to the upper grid value of the cell. The counts in a bin are
the sum over the corners of the counts simulated at the corner times a weight. The
weight is the product over the parameters of ``t_d`` if the corner is at the upper
grid value of parameter ``d``, and of ``1 - t_d`` if it is at the lower one. With one
parameter this is ``(1 - t) n_k + t n_{k+1}``, where ``n_k`` and ``n_{k+1}`` are the
counts at the two grid values around the point.

This works when a parameter changes the counts in each bin smoothly, such as a
dead layer or a source position. It does not work for a parameter that moves a
peak in energy, such as the energy scale: two shifted peaks interpolate to two
peaks of half the height.

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

function GeneralisedHistogram(
    hists::AbstractArray{<:Histogram{<:Real,1},N};
    grid...,
) where {N}
    grid = NamedTuple(grid)
    1 <= N <= 5 ||
        throw(ArgumentError("a generalised histogram has 1 to 5 parameters, not $N"))
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
