"""
    rebin_hist(h::Histogram, factor::Int)

Rebins a `Histogram` by combining every `factor` consecutive bins of a 1D histogram.
If the number of bins is not divisible by `factor`, the last new bin is narrower.
"""
function rebin_hist(h::Histogram{<:Real,1}, factor::Int)
    factor ≥ 1 || throw(ArgumentError("factor must be ≥ 1"))

    weights = h.weights
    edges = h.edges[1]

    new_weights = [sum(weights[i:min(i+factor-1, end)]) for i = 1:factor:length(weights)]

    new_edges = collect(edges[1:factor:end])
    new_edges[end] == edges[end] || push!(new_edges, edges[end])

    return Histogram((new_edges,), new_weights)
end

"""
    rebin_hist(h::Histogram, edges::AbstractVector{<:Real})

Rebins a 1D `Histogram` to the bin edges `edges`, in the units of `h` (keV for a
spectrum), with bins of any width. Each new bin is the sum of the bins of `h` inside it,
so every value in `edges` must be an edge of `h`. The new histogram covers only the range
of `edges`.
"""
function rebin_hist(h::Histogram{<:Real,1}, edges::AbstractVector{<:Real})
    length(edges) >= 2 || throw(ArgumentError("edges needs at least two values"))
    all(>(0), diff(edges)) || throw(ArgumentError("edges must be strictly increasing"))
    weights = [integrate_hist(h, edges[i], edges[i+1]) for i = 1:(length(edges)-1)]
    return Histogram((collect(Float64, edges),), weights)
end

"""
    uniform_binning_with_peaks(h::Histogram, bin_width; peaks = [], blinded, low, high)

Bin edges in keV for the keyword `binning` of [`build_forward_model`](@ref): bins of
`bin_width` keV from `low` to `high`, by default the range of `h`, except that each peak
region `(min, max)` in `peaks` becomes a single bin. A uniform edge inside a peak region,
or less than `bin_width` from it, is dropped, so the continuum bins on both sides of a
peak are between one and two bin widths wide instead of leaving narrow bins. If
`bin_width` does not fit into the range, the last bin is narrower.

Every edge is moved to the nearest edge of `h`, so that the histograms can be rebinned to
them. The peak regions must lie inside the range and must not overlap. For example, with
the ``^{40}``K and ``^{208}``Tl lines each as one bin:

```julia
binning = uniform_binning_with_peaks(h, 20; peaks = [(1455, 1467), (2605, 2625)])
```

`blinded = (min, max)` leaves out a region in keV, for example around the Q-value of
2039 keV. Its edges are treated like those of a peak region, but instead of a bin there
is a gap: the result is then two vectors of edges, below and above the region, which
[`build_forward_model`](@ref) reads as two ranges with the region between them skipped.

The same as `rebin_uniform_with_peaks` of pygama, except that the peaks are given, not
searched for.
"""
function uniform_binning_with_peaks(
    h::Histogram{<:Real,1},
    bin_width::Real;
    peaks = Tuple{Float64,Float64}[],
    blinded::Union{Nothing,Tuple{Real,Real}} = nothing,
    low::Real = first(h.edges[1]),
    high::Real = last(h.edges[1]),
)
    bin_width > 0 || throw(ArgumentError("bin_width must be positive"))
    # the blinded region is placed like a peak region, and its bin is removed at the end
    regions = _binning_regions(peaks, blinded, low, high)

    edges = collect(Float64, low:bin_width:high)
    last(edges) < high && push!(edges, high)
    for (a, b) in regions
        # the uniform edges inside the peak region or less than one bin width from it go,
        # its own two edges come in
        filter!(e -> e <= a - bin_width || e >= b + bin_width || e == a || e == b, edges)
        append!(edges, (a, b))
    end

    fine = h.edges[1]
    edges = unique([_nearest_edge(fine, e) for e in sort!(edges)])
    return _split_at_blinded(edges, fine, blinded)
end

"""
    bayesian_blocks_binning(h::Histogram; p0 = 0.05, peaks = [], blinded, low, high)

Bin edges in keV for the keyword `binning` of [`build_forward_model`](@ref), chosen from
the counts of `h` with the Bayesian blocks method of
[BayesianBlocks.jl](https://github.com/gipert/BayesianBlocks.jl) (Scargle et al. 2013,
<https://doi.org/10.1088/0004-637X/764/2/167>). The bins are as wide as possible while
the counts per keV inside each are consistent with being constant: narrow bins where the
spectrum changes quickly, such as at a peak or a Compton edge, and wide bins on a flat
continuum. Each edge is an edge of `h`. The range is `low` to `high`, by default that of
`h`.

`h` must hold whole numbers of counts, for example the measured spectrum: the method compares
the changes in counts per keV to their Poisson fluctuations. `p0` is roughly the
probability of a false edge, one caused by a fluctuation. A smaller `p0` gives fewer and
wider bins; typical values are 0.01 to 0.1.

Each peak region `(min, max)` in `peaks` becomes a single bin, and the Bayesian blocks are
found separately on each part of the continuum between them, so the counts in a peak do
not change the continuum bins. `blinded = (min, max)` leaves out a region, as in
[`uniform_binning_with_peaks`](@ref): the result is then two vectors of edges, below and
above it, and the counts in the blinded region are not used. All edges are first moved to
the nearest edge of `h`.

```julia
binning = bayesian_blocks_binning(data; peaks = [(1455, 1467), (2605, 2625)])
```

The same as `rebin_bblocks` and `rebin_bblocks_with_peaks` of pygama, except that the
peaks are given, not searched for.
"""
function bayesian_blocks_binning(
    h::Histogram{<:Real,1};
    p0::Real = 0.05,
    peaks = Tuple{Float64,Float64}[],
    blinded::Union{Nothing,Tuple{Real,Real}} = nothing,
    low::Real = first(h.edges[1]),
    high::Real = last(h.edges[1]),
)
    0 < p0 < 1 || throw(ArgumentError("p0 must be between 0 and 1"))
    all(w -> w >= 0 && isinteger(w), h.weights) ||
        throw(ArgumentError("h must hold non-negative whole numbers of counts"))
    regions = _binning_regions(peaks, blinded, low, high)

    # the Bayesian blocks of the continuum from `a` to `b` keV
    fine = h.edges[1]
    function blocks(a, b)
        i, j = _nearest_index(fine, a), _nearest_index(fine, b)
        return _bayesian_blocks(fine[i:j], h.weights[i:(j-1)], p0)
    end
    # each region, a peak or the blinded one, is one bin between two parts of the
    # continuum, each partitioned on its own
    edges = Float64[]
    start = low
    for (a, b) in regions
        append!(edges, blocks(start, a))
        start = b
    end
    append!(edges, blocks(start, high))
    return _split_at_blinded(unique(edges), fine, blinded)
end

# The edges of the Bayesian blocks of the bins with `edges` and `counts`, from
# BayesianBlocks.jl. It places the edges halfway between the centres of bins with counts,
# so each is moved to the nearest edge of the bins, and the outer ones to the ends.
# BayesianBlocks.jl takes the outer blocks to start and end at the centre of the first
# and last bin, half a bin short, so their counts per keV look too high and the first
# and last bin often come out as blocks of their own. Such a block is merged into the
# next one.
function _bayesian_blocks(edges, counts, p0)
    any(>(0), counts) || return [first(edges), last(edges)]
    blocks = BayesianBlocks.partition(Histogram((edges,), round.(Int, counts)); p0)
    inner = unique([_nearest_edge(edges, e) for e in blocks[2:(end-1)]])
    return [first(edges); filter(e -> edges[2] < e < edges[end-1], inner); last(edges)]
end

# the peak regions and the blinded region (min, max) in keV, sorted, checked to lie
# inside (low, high) and not to overlap
function _binning_regions(peaks, blinded, low, high)
    high > low || throw(ArgumentError("high must be above low"))
    regions = [(Float64(a), Float64(b)) for (a, b) in peaks]
    isnothing(blinded) || push!(regions, Float64.(blinded))
    sort!(regions)
    for (a, b) in regions
        low <= a < b <= high ||
            throw(ArgumentError("the region ($a, $b) is not inside ($low, $high)"))
    end
    for (r, next) in zip(regions, regions[2:end])
        r[2] <= next[1] || throw(ArgumentError("the regions $r and $next overlap"))
    end
    return regions
end

# the index of the edge in the sorted `fine` nearest to `e`, the lower one on a tie
function _nearest_index(fine, e)
    i = clamp(searchsortedfirst(fine, e), 2, length(fine))
    return e - fine[i-1] <= fine[i] - e ? i - 1 : i
end

# the edge in the sorted `fine` nearest to `e`
_nearest_edge(fine, e) = fine[_nearest_index(fine, e)]

# `edges` unchanged without a blinded region, otherwise the edges below and above it as
# two ranges, and a part with no bin is left out
function _split_at_blinded(edges, fine, blinded)
    isnothing(blinded) && return edges
    a, b = (_nearest_edge(fine, e) for e in blinded)
    return filter(r -> length(r) >= 2, [filter(<=(a), edges), filter(>=(b), edges)])
end

"""
    add_histograms(hs::Vector{<:Histogram})

Sum the weights of histograms that share the same binning.
"""
function add_histograms(hs::Vector{<:Histogram})
    all(h -> h.edges == hs[1].edges, hs) || error("Histogram binning mismatch")
    return Histogram(hs[1].edges, sum(h.weights for h in hs))
end

"""
    integrate_hist(h::Histogram{<:Any,1}, a, b)

Total counts (sum of weights) of a 1D histogram between `a` and `b`,
which must both lie on bin edges.
"""
function integrate_hist(h::Histogram{<:Any,1}, a, b)
    edges = h.edges[1]

    i_start = searchsortedfirst(edges, a)
    i_end = searchsortedfirst(edges, b)

    (i_start ≤ length(edges) && edges[i_start] ≈ a) ||
        throw(ArgumentError("a = $a must be a bin edge"))
    (i_end ≤ length(edges) && edges[i_end] ≈ b) ||
        throw(ArgumentError("b = $b must be a bin edge"))

    return sum(@view h.weights[i_start:(i_end-1)])
end

"""
    get_valid_bins(edges, range)

Boolean mask of the bins fully inside `range = (min, max)`.
"""
function get_valid_bins(edges::AbstractVector{<:Real}, range::Tuple)
    return (edges[1:(end-1)] .>= range[1]) .& (edges[2:end] .<= range[2])
end

"""
    select_bins(h::Histogram; binning = 1, e_low = -Inf, e_high = Inf)

Rebin the 1D histogram `h` and keep the bins of the fit, in the same way as
[`build_forward_model`](@ref): `binning` rebins it (see [Binning](binning.md)), then only
the bins fully inside `[e_low, e_high]` (keV) are kept, and with several ranges in
`binning` only those inside one of them. By default the whole range of `h` is kept.

Returns `(counts, bin_edges)`: the counts in the kept bins as a flat
`Vector{Float64}`, and `bin_edges = (low = ..., high = ...)`, their lower and upper
edges in keV. The kept bins can have gaps between them, for example a blinded region.
Use it for the measured spectrum and the simulated ones of a forward model written by
hand, with the same keywords, so that they have the same bins.
"""
function select_bins(
    h::Histogram{<:Real,1};
    binning::Union{Int,AbstractVector} = 1,
    e_low::Real = -Inf,
    e_high::Real = Inf,
)
    rebin_to, ranges = _binning(binning)
    rebinned = rebin_hist(h, rebin_to)
    edges = rebinned.edges[1]
    kept = _kept_bins(edges, e_low, e_high, ranges)
    bin_edges = (
        low = collect(Float64, edges[1:(end-1)][kept]),
        high = collect(Float64, edges[2:end][kept]),
    )
    return collect(Float64, rebinned.weights[kept]), bin_edges
end

# true for the bins fully inside [e_low, e_high] and, with ranges, inside one of them
_kept_bins(edges, e_low, e_high, ranges) =
    get_valid_bins(edges, (e_low, e_high)) .& _in_ranges(edges, ranges)

"""
    get_counts(data::Histogram, bin_edges)

Measured counts in the bins of the forward model, as a flat `Vector{Float64}`.

`bin_edges = (low = ..., high = ...)` are the lower and upper edges of the bins, for
example the ones returned by [`build_forward_model`](@ref). The result is not a
histogram: it has one entry per bin, in the same order as `bin_edges`. Each entry is
the sum of the `data` bins between `low[i]` and `high[i]`, which must be edges of
`data`.
"""
function get_counts(data::Histogram{<:Real,1}, bin_edges::NamedTuple{(:low, :high)})
    length(bin_edges.low) == length(bin_edges.high) ||
        throw(DimensionMismatch("bin_edges.low and bin_edges.high differ in length"))
    return Float64[integrate_hist(data, lo, hi) for (lo, hi) in zip(bin_edges...)]
end

"""
    counts_to_histogram(counts, bin_edges)

Turn a flat vector of `counts` back into a `Histogram`, for example to plot the
model prediction. This is the inverse of [`get_counts`](@ref).

`counts[i]` is the content of the bin between `bin_edges.low[i]` and
`bin_edges.high[i]`. The bins must be sorted and must not overlap. The histogram
covers the range from `low[1]` to `high[end]`. A gap between two bins, such as a
region skipped by the binning, becomes one bin with 0 counts.
"""
function counts_to_histogram(
    counts::AbstractVector{<:Real},
    bin_edges::NamedTuple{(:low, :high)},
)
    low, high = bin_edges.low, bin_edges.high
    length(counts) == length(low) == length(high) || throw(
        DimensionMismatch("counts, bin_edges.low and bin_edges.high differ in length"),
    )
    isempty(counts) && throw(ArgumentError("bin_edges is empty"))

    edges = Float64[low[1]]
    weights = Float64[]
    for i in eachindex(counts)
        if i > 1 && !(low[i] ≈ high[i-1])
            low[i] > high[i-1] || throw(ArgumentError("bins are not sorted or overlap"))
            push!(edges, low[i])  # the gap starts at the last edge
            push!(weights, 0.0)
        end
        push!(edges, high[i])
        push!(weights, counts[i])
    end
    return Histogram((edges,), weights)
end
