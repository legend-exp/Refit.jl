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
    uniform_binning_with_peaks(h::Histogram, bin_width; peaks = [], low, high)

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

The same as `rebin_uniform_with_peaks` of pygama, except that the peaks are given, not
searched for.
"""
function uniform_binning_with_peaks(
    h::Histogram{<:Real,1},
    bin_width::Real;
    peaks = Tuple{Float64,Float64}[],
    low::Real = first(h.edges[1]),
    high::Real = last(h.edges[1]),
)
    bin_width > 0 || throw(ArgumentError("bin_width must be positive"))
    high > low || throw(ArgumentError("high must be above low"))
    regions = sort([(Float64(a), Float64(b)) for (a, b) in peaks])
    for (a, b) in regions
        low <= a < b <= high ||
            throw(ArgumentError("the peak region ($a, $b) is not inside ($low, $high)"))
    end
    for (r, next) in zip(regions, regions[2:end])
        r[2] <= next[1] || throw(ArgumentError("the peak regions $r and $next overlap"))
    end

    edges = collect(Float64, low:bin_width:high)
    last(edges) < high && push!(edges, high)
    for (a, b) in regions
        # the uniform edges inside the peak region or less than one bin width from it go,
        # its own two edges come in
        filter!(e -> e <= a - bin_width || e >= b + bin_width || e == a || e == b, edges)
        append!(edges, (a, b))
    end

    # the nearest edge of `h` to each edge
    fine = h.edges[1]
    return unique([fine[argmin(abs.(fine .- e))] for e in sort!(edges)])
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
