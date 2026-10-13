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
    :total in keys(parameters) && throw(
        ArgumentError("a source cannot be called :total, the name is used for the sum"),
    )

    hists = Dict{Symbol,Histogram{Float64,1}}(
        :total => counts_to_histogram(model(parameters), bin_edges),
    )
    for name in keys(parameters)
        # the model with the activities of all the other sources set to zero
        only_one =
            map(k -> k == name ? parameters[k] : zero(parameters[k]), keys(parameters))
        hists[name] =
            counts_to_histogram(model(NamedTuple{keys(parameters)}(only_one)), bin_edges)
    end
    return hists
end
