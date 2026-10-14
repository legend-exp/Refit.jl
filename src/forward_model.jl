"""
    build_forward_model(components; e_low, e_high, blind, norm, binning)

Build a function that maps fit parameters to the predicted counts per energy bin.

`components` is a dictionary `Dict(name => component)` with one simulated 1D energy
spectrum in keV per component. Each component is one of:
- a `Histogram`: its activity is the fit parameter with the same name, `params[name]`
- `(pdf = histogram, activity = f)`: its activity is `f(params)`, a function of the fit
  parameters. For example to fit a half-life instead of an activity, or to tie the
  activities of two components together

`params` is the named tuple of all fit parameters, in any order. The prediction of a
component is `activity * norm * counts` in every bin, and the model is the sum over the
components. Use the global `norm` for unit conversions and the livetime.

The global `norm` multiplies every component, also those with an `activity` function:
their prediction is `f(params) * norm * counts`. So `f` returns an activity, in the same
units as a fit parameter.

All histograms must share the same binning. Only bins fully inside `[e_low, e_high]` and
outside the `blind = (min, max)` region are kept. `binning` merges that many consecutive
bins first.

Returns `(forward_model, bin_edges)`:
- `bin_edges`: `(low, high)`, the lower and upper edges in keV of the kept bins
- `forward_model(params)`: predicted counts per kept bin
- `forward_model(params, Val(:by_component))`: a named tuple with the predicted counts
  of each component, in the order of the sorted names

The time per call grows with the number of bins times the number of components. An
`activity` function adds about 0.1 microseconds per call.
"""
function build_forward_model(
    components::AbstractDict{Symbol};
    e_low::Real = 500,
    e_high::Real = 4000,
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
    binning::Int = 1,
)
    isempty(components) && throw(ArgumentError("components is empty"))
    names = Tuple(sort!(collect(keys(components))))
    specs = map(k -> _component_spec(components[k]), names)
    hists = map(s -> rebin_hist(s.pdf, binning), specs)

    edges = first(hists).edges[1]
    all(h -> h.edges[1] == edges, hists) || error("Histogram binning mismatch")

    is_valid = get_valid_bins(edges, (e_low, e_high), blind)
    bin_low = edges[1:(end-1)][is_valid]
    bin_high = edges[2:end][is_valid]
    counts = map(h -> h.weights[is_valid] .* norm, hists)

    plain, with_function = _split_components(names, specs, counts)
    # the names live in the type of `mc`, so the activities are sorted into the order of
    # the spectra when the model is compiled, not on every call
    mc = NamedTuple{Tuple(plain.names)}(Tuple(plain.counts))
    # all the spectra as a vector, the plain ones first. Julia unrolls a loop over a
    # tuple only up to about 30 entries, a vector has no such limit
    spectra = [plain.counts; with_function.counts]
    order = (plain.names..., with_function.names...)

    # the activity of every component, in the order of `spectra`
    function activities(params)
        fitted = NamedTuple{keys(mc)}(params)
        a = Vector{Float64}(undef, length(spectra))
        for j = 1:length(fitted)
            a[j] = fitted[j]
        end
        for (k, f) in enumerate(with_function.funcs)
            a[length(fitted)+k] = f(params)
        end
        return a
    end

    # the sum of all components
    function forward_model(params::NamedTuple)
        a = activities(params)
        pred = zeros(length(bin_low))
        for j in eachindex(spectra)
            pred .+= a[j] .* spectra[j]
        end
        return pred
    end
    # each component on its own
    function forward_model(params::NamedTuple, ::Val{:by_component})
        a = activities(params)
        parts = NamedTuple{order}(Tuple(a[j] .* spectra[j] for j in eachindex(spectra)))
        return NamedTuple{names}(parts)
    end
    return forward_model, (low = bin_low, high = bin_high)
end

# a component is a Histogram, or a named tuple `(pdf = histogram, activity = f)`
_component_spec(pdf::Histogram{<:Real,1}) = (pdf = pdf, activity = nothing)
function _component_spec(c::NamedTuple)
    extra = setdiff(keys(c), (:pdf, :activity))
    isempty(extra) || throw(ArgumentError("unknown component keys $extra"))
    haskey(c, :pdf) || throw(ArgumentError("a component needs a pdf"))
    return (pdf = c.pdf, activity = get(c, :activity, nothing))
end

"""
    _split_components(names, specs, counts)

Split the components into two lists, each a named tuple of vectors in the same order:
- `plain`: the components whose activity is a fit parameter, `(names, counts)`
- `with_function`: the components with an `activity` function, `(names, funcs, counts)`
"""
function _split_components(names, specs, counts)
    plain = (names = Symbol[], counts = Vector{Float64}[])
    with_function = (names = Symbol[], funcs = Function[], counts = Vector{Float64}[])
    for (name, spec, c) in zip(names, specs, counts)
        if isnothing(spec.activity)
            push!(plain.names, name)
            push!(plain.counts, c)
        else
            push!(with_function.names, name)
            push!(with_function.funcs, spec.activity)
            push!(with_function.counts, c)
        end
    end
    return plain, with_function
end
