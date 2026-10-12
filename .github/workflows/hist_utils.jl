using StatsBase
import Base: +


"""
    rebin_hist(h::Histogram, factor::Int) -> Histogram

Rebins a 1D `StatsBase.Histogram` by combining `factor` consecutive bins.

- Returns a new `Histogram` with rebinned weights and edges.
- Works even if the number of original bins isn’t divisible by `factor`.
- Ensures the constructor uses the correct tuple-of-edges format.
"""
function rebin_hist(h::Histogram{<:Real,1}, factor::Int)
    @assert factor ≥ 1 "Factor must be ≥ 1"

    old_weights = h.weights
    old_edges = h.edges[1]  # extract 1D edges

    # Combine weights
    new_weights =
        [sum(old_weights[i:min(i+factor-1, end)]) for i = 1:factor:length(old_weights)]

    # Create new edges
    new_edges = old_edges[1:factor:end]
    if new_edges[end] != old_edges[end]
        push!(new_edges, old_edges[end])  # ensure last edge matches original
    end

    # Construct a new StatsBase.Histogram properly
    return Histogram((new_edges,), new_weights)
end

"""
    +(h1::Histogram, h2::Histogram) -> Histogram

Adds two histograms with the same binning.
- Returns a new `Histogram` with combined weights.
- Throws an error if the bin edges do not match.

"""
function +(h1::Histogram, h2::Histogram)
    # construct and return the result
    h1.edges == h2.edges || error("Histogram binning mismatch")

    return Histogram(h1.edges, h1.weights .+ h2.weights)

end



"""
    integrate_hist(h::Histogram{<:Any,1}, a, b)

Integrate a 1D histogram over [a, b], assuming `a` and `b`
are exactly bin edges.
"""
function integrate_hist(h::Histogram{<:Any,1}, a, b)
    edges = h.edges[1]
    weights = h.weights

    # Find bin indices
    i_start = searchsortedfirst(edges, a)
    i_end = searchsortedfirst(edges, b) - 1

    # Sanity check (optional but useful)
    @assert edges[i_start] == a "a must align with a bin edge"
    @assert edges[i_end+1] == b "b must align with a bin edge"

    widths = diff(edges)

    total = 0.0
    for i = i_start:i_end
        total += weights[i] * widths[i]
    end

    return total
end

"""
    get_valid_bins(edges::AbstractVector{Real}, range::Tuple, blind::Bool = nothing)

Get bins that are fully inside a specified range and optionally outside a blind region.
- `edges`: Vector of bin edges.
- `range`: Tuple specifying the main range (min, max).
- `blind`: Optional tuple specifying a blind region (min, max) to exclude.
- Returns a boolean array indicating valid bins.
"""
function get_valid_bins(edges::AbstractVector{Real}, range::Tuple, blind::Bool = nothing)
    # Bin edges
    bin_starts = edges[1:(end-1)]
    bin_ends = edges[2:end]

    # Fully inside the main range
    inside_main = (bin_starts .>= range[1]) .& (bin_ends .<= range[2])

    if blind == nothing
        return inside_main
    end
    # Fully outside the central exclusion (optional)
    outside_central = (bin_ends .<= blind[1]) .| (bin_starts .>= blind[2])
    return inside_main .& outside_central

end
