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
    get_valid_bins(edges, range, blind = nothing)

Boolean mask of bins fully inside `range = (min, max)` and, if `blind = (min, max)`
is given, fully outside the blinded region.
"""
function get_valid_bins(edges::AbstractVector{<:Real}, range::Tuple, blind = nothing)
    bin_starts = edges[1:(end-1)]
    bin_ends = edges[2:end]

    inside_main = (bin_starts .>= range[1]) .& (bin_ends .<= range[2])
    isnothing(blind) && return inside_main

    outside_blind = (bin_ends .<= blind[1]) .| (bin_starts .>= blind[2])
    return inside_main .& outside_blind
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
