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

"""
    get_priors(pdfs, data; ref_ranges, e_low, e_high, blind, norm,
               range_factor, peak_factor)

Estimate a range `0 .. a_max` for the activity of each simulated spectrum in
`pdfs` (a `Dict(name => component)`, see [`build_forward_model`](@ref)), based on the
measured spectrum `data`. Components with an `activity` or an `analytic` function are
skipped: the parameters of the function need priors by hand.
Returns `Dict(name => Interval)`. The entries can be replaced by any 1D distribution
from Distributions.jl before calling [`build_prior_dist`](@ref).

These are not real informative priors. They come from the data and only limit the range of the
fit to improve convergence. We can also use external knowledge for informative priors.

`a_max` is the largest activity for which the prediction stays below the data, in
a window: a safety factor times [`get_max_norm`](@ref) divided by `norm`. It is the
smaller of two values:
- the value in the fit range `[e_low, e_high]`, times `range_factor` (default 1.1),
- the value in the reference range of the source, times `peak_factor` (default 1.5).

`ref_ranges` is a `Dict(key => (low, high))` with energy ranges in keV, for example
around a peak. A spectrum uses the range of the longest key that its name starts
with, for example `:K40_hpge_bulk` uses `:K40`. Spectra without a range only use
the fit range. `e_low`, `e_high`, `blind` and `norm` are as in
[`build_forward_model`](@ref).
"""
function get_priors(
    pdfs::AbstractDict{Symbol},
    data::Histogram{<:Real,1};
    ref_ranges::AbstractDict{Symbol,<:Tuple{Real,Real}} = Dict{Symbol,Tuple{Real,Real}}(),
    e_low::Real = 500,
    e_high::Real = 4000,
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
    range_factor::Real = 1.1,
    peak_factor::Real = 1.5,
)
    priors = Dict{Symbol,Any}()
    for (name, component) in pdfs
        spec = _component_spec(component)
        if !isnothing(spec.activity) || !isnothing(spec.analytic)
            kind = isnothing(spec.activity) ? "an analytic" : "an activity"
            @info "get_priors: $name has $kind function, its parameters need priors by hand"
            continue
        end
        pdf = spec.pdf
        a_max = range_factor * get_max_norm(e_low, e_high, data, pdf; blind) / norm
        isfinite(a_max) || throw(ArgumentError("$name has no counts in the fit range"))

        matching =
            filter(k -> startswith(string(name), string(k)), collect(keys(ref_ranges)))
        if !isempty(matching)
            low, high = ref_ranges[argmax(k -> length(string(k)), matching)]
            a_ref = peak_factor * get_max_norm(low, high, data, pdf; blind) / norm
            a_max = min(a_max, a_ref)
        end
        priors[name] = 0.0 .. a_max
    end
    return priors
end

"""
    build_prior_dist(ranges)

Build a prior distribution from the entries of [`get_priors`](@ref), a
`Dict(name => prior)`. Each prior is either a range (`Interval`), which becomes a
uniform distribution, or any 1D distribution from Distributions.jl, which is used as is. It is a `NamedTupleDist` (from ValueShapes.jl) and `rand(prior)` draws a
named tuple of activities, which can be passed to the forward model.
"""
function build_prior_dist(ranges::AbstractDict{Symbol})
    return NamedTupleDist(; (k => _to_distribution(r) for (k, r) in ranges)...)
end

_to_distribution(r::Interval) = Uniform(r.left, r.right)
_to_distribution(d::UnivariateDistribution) = d
