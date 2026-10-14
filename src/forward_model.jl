"""
    build_forward_model(components; e_low, e_high, blind, norm, binning, edges)

Build a function that maps fit parameters to the predicted counts per energy bin.

`components` is a dictionary `Dict(name => component)`. Each component is one of:
- a `Histogram`, a simulated 1D energy spectrum in keV: its activity is the fit
  parameter with the same name, `params[name]`
- `(pdf = histogram, activity = f)`: its activity is `f(params)`, a function of the fit
  parameters. For example to fit a half-life instead of an activity, or to tie the
  activities of two components together
- `(pdf = histogram, analytic = g)`, also with `activity = f`: the simulated spectrum is
  multiplied in every bin by `g(E, params)`, a function of the energy `E` in keV. For
  example an efficiency curve, or the factor ``(Q - E)^2`` of a beta decay
- `(analytic = g,)`: an analytic spectrum with no simulation, for example a flat
  background. `g(E, params)` is the expected number of counts per keV at the energy `E`

`params` is the named tuple of all fit parameters, in any order. The prediction of a
component in bin `i` is `activity * norm * g(E_i, params) * counts[i]`, with `g = 1`
without `analytic`, and `g(E_i, params) * width_i` for an analytic component with no
simulation. `E_i` is the energy at the centre of the bin and `width_i` its width in keV.
The model is the sum over the components. Use the global `norm` for unit conversions and
the livetime: it multiplies every component with a histogram, so an `activity` function
returns an activity in the same units as a fit parameter. It does not multiply an
analytic component with no simulation. `g` is evaluated at the bin centre only, which is
inaccurate for features narrower than about one bin, such as peaks.

All histograms must share the same binning. Without any (only analytic components),
give the bin edges in keV with `edges`. Only bins fully inside `[e_low, e_high]` and
outside the `blind = (min, max)` region are kept. `binning` merges that many consecutive
bins first.

Returns `(forward_model, bin_edges)`:
- `bin_edges`: `(low, high)`, the lower and upper edges in keV of the kept bins
- `forward_model(params)`: predicted counts per kept bin
- `forward_model(params, Val(:by_component))`: a named tuple with the predicted counts
  of each component, in the order of the sorted names

The time per call grows with the number of bins times the number of components. An
`activity` function adds about 0.1 microseconds per call, an `analytic` function a call
of the function in every bin.
"""
function build_forward_model(
    components::AbstractDict{Symbol};
    e_low::Real = 500,
    e_high::Real = 4000,
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
    binning::Int = 1,
    edges::Union{AbstractVector{<:Real},Nothing} = nothing,
)
    isempty(components) && throw(ArgumentError("components is empty"))
    names = Tuple(sort!(collect(keys(components))))
    specs = [_component_spec(components[k]) for k in names]
    # the simulated spectra, rebinned. An analytic component with no simulation has none
    hists = [isnothing(s.pdf) ? nothing : rebin_hist(s.pdf, binning) for s in specs]
    model_edges = _model_edges(hists, edges, binning)

    is_valid = get_valid_bins(model_edges, (e_low, e_high), blind)
    bin_low = model_edges[1:(end-1)][is_valid]
    bin_high = model_edges[2:end][is_valid]
    centre = (bin_low .+ bin_high) ./ 2
    width = bin_high .- bin_low

    # counts per kept bin of every component, in the order of `names`: the simulated
    # spectrum times `norm`, or the bin widths for an analytic component with no
    # simulation, whose function is a number of counts per keV
    counts =
        Vector{Float64}[isnothing(h) ? width : h.weights[is_valid] .* norm for h in hists]

    # the positions in `names` of the components whose activity is a fit parameter, of
    # those with an activity function, and of those with an analytic function
    plain = findall(s -> !isnothing(s.pdf) && isnothing(s.activity), specs)
    with_activity = findall(s -> !isnothing(s.activity), specs)
    with_analytic = findall(s -> !isnothing(s.analytic), specs)
    no_analytic = findall(s -> isnothing(s.analytic), specs)
    activity_funcs = Function[specs[k].activity for k in with_activity]
    analytic_funcs = Function[specs[k].analytic for k in with_analytic]

    # the names live in the type of `mc`, so the fit parameters are sorted into the order
    # of `plain` when the model is compiled, not on every call
    mc = NamedTuple{names[plain]}(Tuple(counts[plain]))

    # the activity of every component: its fit parameter, the value of its activity
    # function, or 1 for an analytic component with no simulation
    function activities(params)
        fitted = NamedTuple{keys(mc)}(params)
        a = ones(length(names))
        for (j, k) in enumerate(plain)
            a[k] = fitted[j]
        end
        for (f, k) in zip(activity_funcs, with_activity)
            a[k] = f(params)
        end
        return a
    end

    # add `a * g(E, params) * c` in every bin, with `E` the bin centre in keV. This is a
    # function of its own so that Julia looks up which function `g` is once, when it is
    # called, and not again in every bin
    function add_analytic!(pred, a, g, c, params)
        for i in eachindex(pred)
            pred[i] += a * g(centre[i], params) * c[i]
        end
        return pred
    end

    # the sum of all components
    function forward_model(params::NamedTuple)
        a = activities(params)
        pred = zeros(length(bin_low))
        for k in no_analytic
            pred .+= a[k] .* counts[k]
        end
        for (g, k) in zip(analytic_funcs, with_analytic)
            add_analytic!(pred, a[k], g, counts[k], params)
        end
        return pred
    end
    # each component on its own
    function forward_model(params::NamedTuple, ::Val{:by_component})
        a = activities(params)
        parts = Vector{Vector{Float64}}(undef, length(names))
        for k in no_analytic
            parts[k] = a[k] .* counts[k]
        end
        for (g, k) in zip(analytic_funcs, with_analytic)
            parts[k] = add_analytic!(zeros(length(bin_low)), a[k], g, counts[k], params)
        end
        return NamedTuple{names}(Tuple(parts))
    end
    return forward_model, (low = bin_low, high = bin_high)
end

# a component is a Histogram, or a named tuple with a `pdf` and/or an `analytic`
# function, and an `activity` function if it has a `pdf`
_component_spec(pdf::Histogram{<:Real,1}) =
    (pdf = pdf, activity = nothing, analytic = nothing)
function _component_spec(c::NamedTuple)
    extra = setdiff(keys(c), (:pdf, :activity, :analytic))
    isempty(extra) || throw(ArgumentError("unknown component keys $extra"))
    spec = (
        pdf = get(c, :pdf, nothing),
        activity = get(c, :activity, nothing),
        analytic = get(c, :analytic, nothing),
    )
    if isnothing(spec.pdf)
        isnothing(spec.analytic) &&
            throw(ArgumentError("a component needs a pdf or an analytic function"))
        isnothing(spec.activity) ||
            throw(ArgumentError("an analytic component with no pdf has no activity"))
    end
    return spec
end

# Bin edges of the model in keV. All the simulated spectra must share them. `edges` is
# treated like a spectrum without counts: it is rebinned and checked in the same way.
function _model_edges(hists, edges, binning)
    binnings = AbstractVector[h.edges[1] for h in hists if !isnothing(h)]
    isnothing(edges) || push!(binnings, rebin_hist(Histogram(edges), binning).edges[1])
    isempty(binnings) && throw(
        ArgumentError("a model with only analytic components needs the bin edges `edges`"),
    )
    allequal(binnings) || error("Histogram binning mismatch")
    return first(binnings)
end
