"""
    get_max_norm(low, high, data, pdf; blind = nothing)

Largest factor by which the simulated spectrum `pdf` can be scaled without exceeding
the measured spectrum `data`: the ratio of the counts of `data` and `pdf` between
`low` and `high` (keV). Only bins fully inside that range and outside the `blind`
range `(min, max)` are used. `data` and `pdf` must have the same binning.

Returns `Inf` if `pdf` has no counts in the range.
"""
function get_max_norm(
    low::Real,
    high::Real,
    data::Histogram{<:Real,1},
    pdf::Histogram{<:Real,1};
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
)
    data.edges[1] == pdf.edges[1] || error("Histogram binning mismatch")
    is_valid = get_valid_bins(data.edges[1], (low, high), blind)
    n_pdf = sum(pdf.weights[is_valid])
    n_pdf > 0 || return Inf
    return sum(data.weights[is_valid]) / n_pdf
end

# Safety factors on the largest activity allowed by the data
const RANGE_FACTOR = 1.1  # in the fit range
const PEAK_FACTOR = 1.5  # in a reference range

"""
    get_priors(pdfs, data; ref_ranges, e_low, e_high, blind, norm)

Guess a range `0 .. a_max` for the activity of each simulated spectrum in
`pdfs` (a `Dict(name => Histogram)`), based on the measured spectrum `data`.
Returns `Dict(name => Interval)`.

These are not real priors. They come from the data and only limit the range of the
fit to improve convergence. Use external knowledge for informative priors.

`a_max` is the largest activity for which the prediction stays below the data, in
a window: the safety factor times [`get_max_norm`](@ref) divided by `norm`. It is the
smaller of two values:
- the value in the fit range `[e_low, e_high]`, with a safety factor of 1.1,
- the value in the reference range of the source, with a safety factor of 1.5.

`ref_ranges` is a `Dict(key => (low, high))` with energy ranges in keV, for example
around a peak. A spectrum uses the range of the longest key that its name starts
with, for example `:K40_hpge_bulk` uses `:K40`. Spectra without a range only use
the fit range. `e_low`, `e_high`, `blind` and `norm` are as in
[`build_forward_model`](@ref).
"""
function get_priors(
    pdfs::AbstractDict{Symbol,<:Histogram{<:Real,1}},
    data::Histogram{<:Real,1};
    ref_ranges::AbstractDict{Symbol,<:Tuple{Real,Real}} = Dict{Symbol,Tuple{Real,Real}}(),
    e_low::Real = 500,
    e_high::Real = 4000,
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
)
    priors = Dict{Symbol,Interval{:closed,:closed,Float64}}()
    for (name, pdf) in pdfs
        a_max = RANGE_FACTOR * get_max_norm(e_low, e_high, data, pdf; blind) / norm
        isfinite(a_max) || throw(ArgumentError("$name has no counts in the fit range"))

        matching =
            filter(k -> startswith(string(name), string(k)), collect(keys(ref_ranges)))
        if !isempty(matching)
            low, high = ref_ranges[argmax(k -> length(string(k)), matching)]
            a_ref = PEAK_FACTOR * get_max_norm(low, high, data, pdf; blind) / norm
            a_max = min(a_max, a_ref)
        end
        priors[name] = 0.0 .. a_max
    end
    return priors
end
