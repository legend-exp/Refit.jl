"""
    get_max_norm(low, high, data, pdf; blind = nothing, norm = 1, factor = 1)

Largest activity of the simulated spectrum `pdf` that the `data` allow: the activity
for which the prediction has `factor` times the counts of the data between `low` and
`high` (keV). Only bins fully inside that range and outside the `blind` range
`(min, max)` are used. `norm` is the factor multiplying the spectrum in the forward
model, see [`build_forward_model`](@ref). `data` and `pdf` must have the same binning.

Returns `Inf` if `pdf` has no counts in the range.
"""
function get_max_norm(
    low::Real,
    high::Real,
    data::Histogram{<:Real,1},
    pdf::Histogram{<:Real,1};
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
    factor::Real = 1,
)
    data.edges[1] == pdf.edges[1] || error("Histogram binning mismatch")
    is_valid = get_valid_bins(data.edges[1], (low, high), blind)
    n_pdf = sum(pdf.weights[is_valid])
    n_pdf > 0 || return Inf
    return factor * sum(data.weights[is_valid]) / (n_pdf * norm)
end

"""
    get_priors(pdfs, data; ref_regions, e_low, e_high, blind, norm,
               peak_width, peak_factor, range_factor)

Guess a range `0 .. a_max` for the activity of each simulated spectrum in
`pdfs` (a `Dict(name => Histogram)`), based on the measured spectrum `data`.
Returns `Dict(name => Interval)`.

These are not real priors. They come from the data and only limit the range of the
fit to improve convergence. Use external knowledge for informative priors.

`a_max` is the largest activity for which the prediction stays below the data
(see [`get_max_norm`](@ref)) in two windows. The smaller value is used:
- the fit range `[e_low, e_high]`, with the factor `range_factor` (default 1.1),
- the reference region of the source, with the factor `peak_factor` (default 1.5).

`ref_regions` is a `Dict(key => region)`. A `Real` region is the energy in keV of a
peak, the window is the peak `± peak_width` (default 5 keV). A `(low, high)` tuple
is used as the window. A spectrum uses the region of the longest key that its name
starts with, for example `:K40_hpge_bulk` uses `:K40`. Spectra without a region
only use the fit range. `e_low`, `e_high`, `blind` and `norm` are as in
[`build_forward_model`](@ref).
"""
function get_priors(
    pdfs::AbstractDict{Symbol,<:Histogram{<:Real,1}},
    data::Histogram{<:Real,1};
    ref_regions::AbstractDict{Symbol} = Dict{Symbol,Real}(),
    e_low::Real = 500,
    e_high::Real = 4000,
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
    peak_width::Real = 5,
    peak_factor::Real = 1.5,
    range_factor::Real = 1.1,
)
    priors = Dict{Symbol,Interval{:closed,:closed,Float64}}()
    for (name, pdf) in pdfs
        a_max = get_max_norm(e_low, e_high, data, pdf; blind, norm, factor = range_factor)
        isfinite(a_max) || throw(ArgumentError("$name has no counts in the fit range"))

        region = find_region(name, ref_regions)
        if !isnothing(region)
            window = region isa Real ? (region - peak_width, region + peak_width) : region
            a_peak = get_max_norm(window..., data, pdf; blind, norm, factor = peak_factor)
            a_max = min(a_max, a_peak)
        end
        priors[name] = 0.0 .. a_max
    end
    return priors
end

# Region of the longest key of `ref_regions` that `name` starts with, or `nothing`.
function find_region(name::Symbol, ref_regions::AbstractDict{Symbol})
    keys_matching =
        filter(k -> startswith(string(name), string(k)), collect(keys(ref_regions)))
    isempty(keys_matching) && return nothing
    region = ref_regions[argmax(k -> length(string(k)), keys_matching)]
    (region isa Real || region isa Tuple{Real,Real}) ||
        throw(ArgumentError("region of $name must be a number or a (low, high) tuple"))
    return region
end
