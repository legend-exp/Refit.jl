using LegendHDF5IO
using StatsBase
using YAML

"""
    get_best_fit_histograms(parameters, model, bin_edges)

Histograms of the prediction of the model for the `parameters` (a named tuple of
activities, for example the best fit). `model` and `bin_edges` come from
[`build_forward_model`](@ref).

Returns a `Dict(name => Histogram)` with the total prediction in `:total` and the
contribution of each source under the name of its activity. The blinded region and
the range outside of the fit are empty bins (see [`counts_to_histogram`](@ref)).
"""
function get_best_fit_histograms(
    parameters::NamedTuple,
    model::Function,
    bin_edges::NamedTuple{(:low, :high)},
)
    :total in keys(parameters) &&
        throw(ArgumentError("a source cannot be called :total, the name is used for the sum"))

    hists = Dict{Symbol,Histogram{Float64,1}}(
        :total => counts_to_histogram(model(parameters), bin_edges),
    )
    for name in keys(parameters)
        # the model with the activities of all the other sources set to zero
        only_one = map(k -> k == name ? parameters[k] : zero(parameters[k]), keys(parameters))
        hists[name] =
            counts_to_histogram(model(NamedTuple{keys(parameters)}(only_one)), bin_edges)
    end
    return hists
end

"""
    write_best_fit(file, hists; data = nothing)

Write the histograms of [`get_best_fit_histograms`](@ref) to the lh5 file `file`: the
total in `best_fit/total` and each source in `best_fit/<name>`. An existing file is
overwritten.

With `data` (a `Histogram`, for example from [`counts_to_histogram`](@ref)) the
measured spectrum is also written, in `data`. The histograms can be read back with
[`read_hist`](@ref).
"""
function write_best_fit(
    file::AbstractString,
    hists::AbstractDict{Symbol,<:Histogram{<:Real,1}};
    data::Union{Nothing,Histogram{<:Real,1}} = nothing,
)
    lh5open(file, "w") do f
        for (name, h) in hists
            f["best_fit/$name"] = h
        end
        isnothing(data) || (f["data"] = data)
    end
    return file
end

"""
    get_intervals(samples; level = 0.68)

Central credible interval of each parameter from the posterior `samples`, as a
`Dict(name => (low, high))`. `samples` is the result of `bat_sample`, with the
parameters in `samples.v` and the weights in `samples.weight`. The interval contains
the fraction `level` of the posterior, with the same probability below and above it.
For `level = 0.68`, `low` and `high` are the 16% and 84% quantiles.
"""
function get_intervals(samples; level::Real = 0.68)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    weights = Weights(samples.weight)
    probabilities = [(1 - level) / 2, (1 + level) / 2]
    return Dict(
        name => Tuple(quantile([v[name] for v in samples.v], weights, probabilities)) for
        name in keys(first(samples.v))
    )
end

"""
    write_activities(file, parameters, samples; level = 0.68)

Write the activities to the YAML file `file`: for each one the best fit value (from
`parameters`, a named tuple, for example `mode(samples)`) and the central credible
interval `low` to `high` from [`get_intervals`](@ref).
"""
function write_activities(
    file::AbstractString,
    parameters::NamedTuple,
    samples;
    level::Real = 0.68,
)
    intervals = get_intervals(samples; level)
    activities = Dict(
        string(name) => Dict(
            "best_fit" => Float64(parameters[name]),
            "low" => intervals[name][1],
            "high" => intervals[name][2],
        ) for name in keys(parameters)
    )
    YAML.write_file(file, activities)
    return file
end
