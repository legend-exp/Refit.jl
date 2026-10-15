"""
    GeneralisedHistogram(hists, grid...)

A histogram whose counts depend on parameters. It holds the histogram simulated at
every point of a grid in the parameters, and interpolates linearly between the grid
points, in each parameter and separately in each bin. [`get_weights`](@ref) and
[`get_histogram`](@ref) return the counts at given values of the parameters.

`hists` is an array of histograms with one dimension per parameter, all with the same
binning. Each argument of `grid` is a pair `name => values` that describes one dimension
of `hists`, in order: the first pair the first dimension, the second pair the second,
and so on. The values are an evenly spaced range of at least two values. For example

```julia
hists = [simulate(fccd, tl) for fccd in 0.6:0.2:1.4, tl in 0.0:0.1:0.5]
ghist = GeneralisedHistogram(hists, :fccd => 0.6:0.2:1.4, :tl => 0.0:0.1:0.5)
```

gives `hists[i, j]` the `i`-th FCCD value and the `j`-th TL value. If two grids have the
same length, only the order of the pairs tells them apart.

Each evaluation sums the histograms of the ``2^N`` corners of the grid cell around the
parameter values, so its time doubles with each parameter.

!!! warning
    The interpolation is an approximation. It suits a parameter that changes the counts
    in each bin smoothly. It does not suit a parameter that moves a peak, such as the
    energy scale: between two grid points the peak appears twice, at half the height.
    Check that the grid is fine enough for each use.

# Fields
* `bin_edges`: `(low = ..., high = ...)`, the lower and upper edges of the bins, in the
  units of the histograms (keV for a spectrum). After [`select_bins`](@ref) the bins
  can have gaps between them
* `grid`: a named tuple with the grid values of each parameter
* `counts`: an array with the counts per bin at each grid point, `counts[i, j]`
  belongs to the `i`-th value of the first parameter and the `j`-th of the second
"""
struct GeneralisedHistogram{N,G<:NamedTuple}
    bin_edges::@NamedTuple{low::Vector{Float64}, high::Vector{Float64}}
    grid::G
    counts::Array{Vector{Float64},N}
end

function GeneralisedHistogram(
    hists::AbstractArray{<:Histogram{<:Real,1},N},
    grid::Pair{Symbol}...,
) where {N}
    length(grid) >= 1 || throw(ArgumentError("a generalised histogram needs a parameter"))
    length(grid) == N || throw(
        ArgumentError("hists has $N dimensions but there are $(length(grid)) parameters"),
    )
    for (name, values) in grid
        values isa AbstractRange ||
            throw(ArgumentError("the grid of $name must be a range, such as 0.5:0.1:2.0"))
        length(values) >= 2 && step(values) > 0 ||
            throw(ArgumentError("the grid of $name needs at least two increasing values"))
    end
    grid = NamedTuple(grid)
    size(hists) == map(length, Tuple(grid)) || throw(
        DimensionMismatch(
            "hists has size $(size(hists)), the grid $(map(length, Tuple(grid)))",
        ),
    )

    edges = collect(Float64, first(hists).edges[1])
    all(h -> h.edges[1] == edges, hists) ||
        throw(ArgumentError("all histograms must have the same binning"))
    counts = map(h -> collect(Float64, h.weights), hists)
    return GeneralisedHistogram(_bin_edges(edges), grid, counts)
end

# the lower and upper edges of the bins of a vector of edges
_bin_edges(edges) =
    (low = collect(Float64, edges[1:(end-1)]), high = collect(Float64, edges[2:end]))

# the edges of bins without gaps between them, or an error
function _contiguous_edges(bin_edges)
    bin_edges.low[2:end] == bin_edges.high[1:(end-1)] || throw(
        ArgumentError("the bins have gaps: rebin before selecting bins with select_bins"),
    )
    return [bin_edges.low; bin_edges.high[end]]
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
    weights = zeros(promote_type(Float64, map(typeof, t)...), length(gh.bin_edges.low))
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
A gap between the bins, left by [`select_bins`](@ref), becomes one bin with 0 counts
(see [`counts_to_histogram`](@ref)).
"""
get_histogram(gh::GeneralisedHistogram, params::NamedTuple) =
    counts_to_histogram(get_weights(gh, params), gh.bin_edges)

"""
    rebin_hist(gh::GeneralisedHistogram, binning)

Rebin the spectrum at every grid point of `gh`, with an integer number of bins to
merge or a vector of new bin edges, as for a `Histogram`. Rebinning and interpolation
are both sums of counts, so rebinning before interpolating gives the same result as
interpolating before rebinning. The bins of `gh` must not have gaps.
"""
function rebin_hist(gh::GeneralisedHistogram, binning::Union{Int,AbstractVector{<:Real}})
    edges = _contiguous_edges(gh.bin_edges)
    hists = map(c -> rebin_hist(Histogram((edges,), c), binning), gh.counts)
    return GeneralisedHistogram(hists, pairs(gh.grid)...)
end

"""
    select_bins(gh::GeneralisedHistogram; binning = 1, e_low = -Inf, e_high = Inf)

The same as [`select_bins`](@ref) for a `Histogram`, at every grid point of `gh`.
Returns `(gh_selected, bin_edges)`: a generalised histogram with only the selected bins,
and their lower and upper edges in keV.
"""
function select_bins(
    gh::GeneralisedHistogram;
    binning::Union{Int,AbstractVector} = 1,
    e_low::Real = -Inf,
    e_high::Real = Inf,
)
    rebin_to, ranges = _binning(binning)
    rebinned = rebin_hist(gh, rebin_to)
    edges = _contiguous_edges(rebinned.bin_edges)
    kept = _kept_bins(edges, e_low, e_high, ranges)
    bin_edges = (low = edges[1:(end-1)][kept], high = edges[2:end][kept])
    selected = GeneralisedHistogram(bin_edges, gh.grid, map(c -> c[kept], rebinned.counts))
    return selected, bin_edges
end

"""
    add_histograms(ghs::Vector{<:GeneralisedHistogram})

Sum the counts of generalised histograms with the same binning and the same grid,
grid point by grid point.
"""
function add_histograms(ghs::Vector{<:GeneralisedHistogram})
    all(gh -> gh.bin_edges == ghs[1].bin_edges, ghs) || error("Histogram binning mismatch")
    all(gh -> gh.grid == ghs[1].grid, ghs) || error("Grid mismatch")
    return GeneralisedHistogram(ghs[1].bin_edges, ghs[1].grid, sum(gh.counts for gh in ghs))
end
