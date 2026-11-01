"""
    activity_prior(; val, unc, unit = u"μBq")
    activity_prior(; limit, cl, unit = u"μBq")

Prior distribution of an activity, from a measured value `val` with its uncertainty `unc`
(one standard deviation), or from an upper limit `limit` at the confidence level `cl`,
for example 0.9. `val`, `unc` and `limit` are activities with units (from Unitful.jl),
for example `83u"μBq"`.

A measured value gives a Gaussian cut at 0, `truncated(Normal(val, unc); lower = 0)`. An
upper limit gives a half-normal distribution: a Gaussian centred at 0 and cut at 0,
whose quantile `cl` is `limit`. The distribution is of the activity in `unit`, as a
number: the sources of [`build_forward_model`](@ref) need numbers, and their `norm` must
turn this unit into the number of decays, for example `ustrip(NoUnits, unit * livetime)`.
"""
function activity_prior(;
    val::Union{Nothing,Unitful.Quantity} = nothing,
    unc::Union{Nothing,Unitful.Quantity} = nothing,
    limit::Union{Nothing,Unitful.Quantity} = nothing,
    cl::Union{Nothing,Real} = nothing,
    unit::Unitful.Units = u"μBq",
)
    given = .!isnothing.((val, unc, limit, cl))
    μ, σ = if given == (true, true, false, false)
        ustrip(unit, val), ustrip(unit, unc)
    elseif given == (false, false, true, true)
        0 < cl < 1 || throw(ArgumentError("cl must be between 0 and 1"))
        # a Gaussian of width σ centred at 0 is within ±limit with probability cl
        0, ustrip(unit, limit) / quantile(Normal(), (1 + cl) / 2)
    else
        throw(ArgumentError("give either val and unc, or limit and cl"))
    end
    return truncated(Normal(μ, σ); lower = 0)
end

"""
    read_activities(file; unit = u"μBq")
    read_activities(prod::SimProd; unit = u"μBq")

Read the activities of the components of an experiment from the YAML file `file`, or from
`bkgmodel/<experiment>/components.yaml` in the configuration directory of the production
cycle `prod`.

The file maps the name of each component to its quantity and to the activity of each
simulation (`simid`) that stands for it. The quantity is one of `mass_in_kg`,
`area_in_m2` or `number_of_pieces`, with the activity per unit of quantity in
`specific_activity_in_uBq_per_kg`, `specific_activity_in_uBq_per_m2` or
`activity_per_piece_in_uBq`, in the units of their names. Each activity is
`{ val, unc }` or `{ limit, cl }`, as in [`activity_prior`](@ref):

```yaml
wls_reflector:
  area_in_m2: 36
  specific_activity_in_uBq_per_m2:
    wls_reflector_Pb214_to_Po214: { val: 78, unc: 30 }
    wls_reflector_Pb212_to_Pb208: { limit: 12, cl: 0.9 }
```

Returns `Dict(component => Dict(simid => distribution))`, with the names as `Symbol`s and
the distribution of the activity of the component in `unit` (see
[`activity_prior`](@ref)): the quantity times the activity per unit.
"""
function read_activities(file::AbstractString; unit::Unitful.Units = u"μBq")
    # the field of the quantity and of the activity per unit of quantity, whose product is
    # an activity in µBq
    fields = Dict(
        "mass_in_kg" => "specific_activity_in_uBq_per_kg",
        "area_in_m2" => "specific_activity_in_uBq_per_m2",
        "number_of_pieces" => "activity_per_piece_in_uBq",
    )
    activities = Dict{Symbol,Dict{Symbol,UnivariateDistribution}}()
    for (component, entry) in YAML.load_file(file)
        found = filter(q -> haskey(entry, q), collect(keys(fields)))
        length(found) == 1 && haskey(entry, fields[only(found)]) ||
            throw(ArgumentError("$component in $file needs one quantity with its activity"))
        q = entry[only(found)]
        activities[Symbol(component)] = Dict(
            Symbol(simid) => activity_prior(;
                (Symbol(k) => (k == "cl" ? v : q * v * u"μBq") for (k, v) in a)...,
                unit,
            ) for (simid, a) in entry[fields[only(found)]]
        )
    end
    return activities
end

function read_activities(prod::SimProd; unit::Unitful.Units = u"μBq")
    config = YAML.load_file(joinpath(prod.path, "simflow-config.yaml"))
    config_dir = _cycle_path(prod.path, config["paths"]["config"])
    return read_activities(
        joinpath(config_dir, "bkgmodel", prod.experiment, "components.yaml");
        unit,
    )
end

"""
    get_activity(activities, simid)

The prior of the activity of the simulation `simid` in the `activities` of
[`read_activities`](@ref), for example to replace the range of [`get_priors`](@ref) for
the parameter of its spectrum. It is an error if no component lists `simid`, or if
several do: the activity of the simid is then their sum, which is left to the caller.
"""
function get_activity(activities::AbstractDict{Symbol}, simid::Symbol)
    found = [component for (component, a) in activities if haskey(a, simid)]
    isempty(found) && throw(ArgumentError("no component has an activity for $simid"))
    length(found) == 1 || throw(
        ArgumentError(
            "several components have an activity for $simid: $(join(found, ", "))",
        ),
    )
    return activities[only(found)][simid]
end
