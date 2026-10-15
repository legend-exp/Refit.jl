"""
    build_forward_model(components; e_low, e_high, norm, binning, edges, points_per_bin)

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

In place of a `Histogram`, a component can have a [`GeneralisedHistogram`](@ref), a
simulated spectrum that depends on fit parameters, for example the thickness of a dead
layer. At each evaluation it is interpolated at the values of these parameters in
`params`, then treated like a `Histogram`.

`params` is the named tuple of all fit parameters, in any order. The prediction of a
component in bin `i` is `activity * norm * g(E_i, params) * counts[i]`, with `g = 1`
without `analytic`, and `g(E_i, params) * width_i` for an analytic component with no
simulation. `E_i` is the energy at the centre of the bin and `width_i` its width in keV.
The model is the sum over the components. Use the global `norm` for unit conversions and
the livetime: it multiplies every component with a histogram, so an `activity` function
returns an activity in the same units as a fit parameter. It does not multiply an
analytic component with no simulation. By default `g` is evaluated at the bin centre
only, which is inaccurate for features narrower than about one bin, such as peaks. With
`points_per_bin = n`, the bin is split into `n` equal parts, `g` is evaluated at the
centre of each and the `n` values are averaged. This takes `n` times as many calls of
`g`, and the error of the integral over the bin falls as ``1/n^2``.

All histograms must share the same binning. Without any (only analytic components),
give the bin edges in keV with `edges`. The histograms are first rebinned with
`binning`, for example read from a YAML file:
- an integer merges that many consecutive bins
- a vector of numbers gives the new bin edges in keV, with bins of any width
- a vector of ranges gives each range its own bin edges in keV, for example
  `[100:20:3000, 3200:100:4200, [4200, 4250, 4270]]`. A range is a vector or range of
  edges, or a string `"start:step:stop"` as in a YAML file. The regions between the
  ranges, here 3000 to 3200 keV, are skipped. Ranges that touch, here at 4200 keV, are
  joined. A string always ends at `stop`: `"500:200:2015"` ends with a narrower bin from
  1900 to 2015 keV. A Julia range stops at its last step instead: `500:200:2015` is
  `500:200:1900`

Each new edge must be an edge of the histograms. Then only the bins fully inside
`[e_low, e_high]` are kept. To blind a region, leave it out between two ranges: for
example `binning = [100:20:2000, 2080:20:4000]` does not fit 2000 to 2080 keV, around the
Q-value of 2039 keV.

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
    norm::Real = 1,
    binning::Union{Int,AbstractVector} = 1,
    edges::Union{AbstractVector{<:Real},Nothing} = nothing,
    points_per_bin::Int = 1,
)
    isempty(components) && throw(ArgumentError("components is empty"))
    points_per_bin >= 1 || throw(ArgumentError("points_per_bin must be at least 1"))
    names = Tuple(sort!(collect(keys(components))))
    specs = [_component_spec(components[k]) for k in names]
    for (name, s) in zip(names, specs)
        slow = isnothing(s.analytic) ? Symbol[] : _slow_globals(s.analytic)
        isempty(slow) || @warn "The analytic function of $name reads the global \
            variable(s) $(join(slow, ", ")), neither `const` nor of a fixed type. This \
            makes every call of the model much slower: declare them `const`, or with a \
            type such as `x::Float64 = 1.0`"
    end
    # the counts of the simulated spectra in the bins of the fit, and the edges of these
    # bins. An analytic component with no simulation has none
    selected = [
        isnothing(s.pdf) ? nothing : select_bins(s.pdf; binning, e_low, e_high) for
        s in specs
    ]
    bin_low, bin_high = _model_bin_edges(selected, edges; binning, e_low, e_high)
    width = bin_high .- bin_low
    # the energies in keV at which analytic functions are evaluated, `points_per_bin` of
    # them in every bin: the centres of that many equal parts of the bin
    energies = [
        bin_low[i] + (j - 0.5) * width[i] / points_per_bin for
        j = 1:points_per_bin, i in eachindex(bin_low)
    ]

    # counts per kept bin of every component, in the order of `names`: the simulated
    # spectrum times `norm`, or the bin widths for an analytic component with no
    # simulation, whose function is a number of counts per keV. A generalised histogram
    # is kept in `ghists` instead, and interpolated at every evaluation
    is_ghist = [!isnothing(sel) && first(sel) isa GeneralisedHistogram for sel in selected]
    counts = Vector{Float64}[
        isnothing(sel) ? width : is_ghist[k] ? Float64[] : first(sel) .* norm for
        (k, sel) in enumerate(selected)
    ]
    ghists = Dict{Int,GeneralisedHistogram}(
        k => _times(first(selected[k]), norm) for k in findall(is_ghist)
    )

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

    # add `a * g * c` in every bin, with `g` the average of `g(E, params)` over the
    # energies of the bin. This is a function of its own so that Julia looks up which
    # function `g` is once, when it is called, and not again in every bin
    function add_analytic!(pred, a, g, c, params)
        for i in eachindex(pred)
            g_sum = 0.0
            for j = 1:points_per_bin
                g_sum += g(energies[j, i], params)
            end
            pred[i] += a * g_sum / points_per_bin * c[i]
        end
        return pred
    end

    # the sum of all components
    function forward_model(params::NamedTuple)
        a = activities(params)
        pred = zeros(length(bin_low))
        for k in no_analytic
            if is_ghist[k]
                _add_interpolated!(pred, a[k], ghists[k], params)
            else
                pred .+= a[k] .* counts[k]
            end
        end
        for (g, k) in zip(analytic_funcs, with_analytic)
            if is_ghist[k]
                add_analytic!(pred, a[k], g, get_weights(ghists[k], params), params)
            else
                add_analytic!(pred, a[k], g, counts[k], params)
            end
        end
        return pred
    end
    # each component on its own
    function forward_model(params::NamedTuple, ::Val{:by_component})
        a = activities(params)
        parts = Vector{Vector{Float64}}(undef, length(names))
        for k in no_analytic
            c = is_ghist[k] ? get_weights(ghists[k], params) : counts[k]
            parts[k] = a[k] .* c
        end
        for (g, k) in zip(analytic_funcs, with_analytic)
            c = is_ghist[k] ? get_weights(ghists[k], params) : counts[k]
            parts[k] = add_analytic!(zeros(length(bin_low)), a[k], g, c, params)
        end
        return NamedTuple{names}(Tuple(parts))
    end
    return forward_model, (low = bin_low, high = bin_high)
end

# add `a` times the counts of the generalised histogram `gh` at the parameters. This is
# a function of its own so that Julia knows the type of `gh` inside it
function _add_interpolated!(pred, a, gh::GeneralisedHistogram, params)
    pred .+= a .* get_weights(gh, params)
    return pred
end

# the generalised histogram `gh` with its counts times `norm`
_times(gh::GeneralisedHistogram, norm) =
    GeneralisedHistogram(gh.bin_edges, gh.grid, map(c -> c .* norm, gh.counts))

# a component is a Histogram or a GeneralisedHistogram, or a named tuple with a `pdf`
# and/or an `analytic` function, and an `activity` function if it has a `pdf`
_component_spec(pdf::Union{Histogram{<:Real,1},GeneralisedHistogram}) =
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

# `binning` as what the histograms are rebinned with, an integer or bin edges in keV, and
# the ranges whose bins are kept, `nothing` to keep all. A vector of numbers is a list of
# edges. Otherwise it is a list of ranges, each with its own edges: the histograms are
# rebinned to the edges of all of them, and the bins between two ranges are skipped.
# The vectors may come from a file, so they are recognised by their content, not type.
_binning(binning::Int) = (binning, nothing)
function _binning(binning::AbstractVector)
    all(x -> x isa Real, binning) && return (collect(Float64, binning), nothing)
    ranges = sort([_range_edges(r) for r in binning]; by = first)
    all(r -> length(r) >= 2, ranges) ||
        throw(ArgumentError("each range of `binning` needs at least two edges"))
    for (r, next) in zip(ranges, ranges[2:end])
        last(r) <= first(next) || throw(ArgumentError("the ranges of `binning` overlap"))
    end
    return unique(reduce(vcat, ranges)), ranges
end

# the bin edges in keV of one range of `binning`: a vector or range of numbers, or a
# string "start:step:stop", the way a range is written in a YAML file. A string always
# ends at `stop`: if the steps do not reach it, the last bin is narrower
_range_edges(r::AbstractVector) = collect(Float64, r)
function _range_edges(r::AbstractString)
    parts = split(r, ":")
    length(parts) == 3 ||
        throw(ArgumentError("the range \"$r\" of `binning` is not \"start:step:stop\""))
    start, step, stop = parse.(Float64, parts)
    edges = collect(start:step:stop)
    # `stop` itself, also when rounding left the last step a hair below it
    if last(edges) ≈ stop
        edges[end] = stop
    else
        push!(edges, stop)
    end
    return edges
end
_range_edges(r) =
    throw(ArgumentError("the range $r of `binning` is not a vector or a string"))

# true for the bins inside one of the `ranges`, for all bins without ranges
_in_ranges(edges, ::Nothing) = trues(length(edges) - 1)
_in_ranges(edges, ranges) = [
    any(r -> first(r) <= low && high <= last(r), ranges) for
    (low, high) in zip(edges[1:(end-1)], edges[2:end])
]

# The lower and upper edges in keV of the bins of the model. All the simulated spectra
# must have the same bins. `edges` is treated like a spectrum without counts: its bins
# are selected and checked in the same way.
function _model_bin_edges(selected, edges; binning, e_low, e_high)
    all_bin_edges = @NamedTuple{low::Vector{Float64}, high::Vector{Float64}}[
        last(sel) for sel in selected if !isnothing(sel)
    ]
    if !isnothing(edges)
        empty = Histogram((collect(Float64, edges),), zeros(length(edges) - 1))
        push!(all_bin_edges, last(select_bins(empty; binning, e_low, e_high)))
    end
    isempty(all_bin_edges) && throw(
        ArgumentError("a model with only analytic components needs the bin edges `edges`"),
    )
    allequal(all_bin_edges) || error("Histogram binning mismatch")
    return first(all_bin_edges)
end

# The global variables that `f` reads directly and that are neither `const` nor of a
# fixed type. Julia cannot know their type when it compiles `f`, so every call of `f`
# looks them up and converts them, which is slow. Empty if the code of `f` cannot be
# read, which uses internals of Julia that may change.
function _slow_globals(f)
    try
        refs = GlobalRef[]
        for m in methods(f)
            foreach(x -> _globalrefs!(refs, x), Base.uncompressed_ast(m).code)
        end
        slow = filter(refs) do g
            isdefined(g.mod, g.name) &&
                !isconst(g.mod, g.name) &&
                Core.get_binding_type(g.mod, g.name) === Any
        end
        return unique(map(g -> g.name, slow))
    catch
        return Symbol[]
    end
end
_globalrefs!(refs, x) = refs
_globalrefs!(refs, x::GlobalRef) = push!(refs, x)
_globalrefs!(refs, x::Expr) = (foreach(a -> _globalrefs!(refs, a), x.args); refs)
