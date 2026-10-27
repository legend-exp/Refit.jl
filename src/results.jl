"""
    get_best_fit_histograms(parameters, model, bin_edges)

Histograms of the prediction of the model for the `parameters` (a named tuple, for
example the best fit). `model` and `bin_edges` come from
[`build_forward_model`](@ref).

Returns an ordered dictionary `name => Histogram` with the total prediction in `:total`,
followed by the contribution of each component under its name in the model, sorted by
name. A region skipped by the binning, such as a blinded one, and the range outside of
the fit are empty bins (see [`counts_to_histogram`](@ref)).
"""
function get_best_fit_histograms(
    parameters::NamedTuple,
    model::Function,
    bin_edges::NamedTuple{(:low, :high)},
)
    parts = model(parameters, Val(:by_component))
    :total in keys(parts) && throw(
        ArgumentError("a component cannot be called :total, the name is used for the sum"),
    )

    hists = OrderedDict{Symbol,Histogram{Float64,1}}(
        :total => counts_to_histogram(model(parameters), bin_edges),
    )
    for name in keys(parts)
        hists[name] = counts_to_histogram(parts[name], bin_edges)
    end
    return hists
end

"""
    get_posterior_bands(samples, model, bin_edges; n = 10_000, groups = nothing,
                        quantiles = (0.16, 0.84))

Uncertainty bands of the prediction of the model, from the posterior `samples` (for
example from `bat_sample`: the parameters of each sample are in `samples.v` and its
weight in `samples.weight`). `model` and `bin_edges` come from
[`build_forward_model`](@ref).

`n` samples are drawn from the posterior, with a probability proportional to their
weight, and the prediction of each one is computed as in
[`get_best_fit_histograms`](@ref). With `groups`, the sources are summed into groups as
in [`group_histograms`](@ref). In each bin, the band goes from the lower to the upper of
the `quantiles` of the `n` predicted counts. The default is the central 68% interval.

Returns an ordered dictionary with the same keys as [`get_best_fit_histograms`](@ref) (or
[`group_histograms`](@ref) with `groups`), with `(low, high)` for each key: the two edges
of the band as `Histogram`s. Draw them with the `band` keyword of [`plot_hist!`](@ref).
"""
function get_posterior_bands(
    samples,
    model::Function,
    bin_edges::NamedTuple{(:low, :high)};
    n::Integer = 10_000,
    groups::Union{Nothing,AbstractVector{<:Pair}} = nothing,
    quantiles::Tuple{Real,Real} = (0.16, 0.84),
)
    n > 0 || throw(ArgumentError("n must be positive"))
    0 <= quantiles[1] < quantiles[2] <= 1 ||
        throw(ArgumentError("quantiles must be increasing and between 0 and 1"))

    draws = sample(eachindex(samples.v), Weights(Float64.(samples.weight)), n)

    # the predicted counts of each source: one row per drawn sample, one column per bin
    predictions = OrderedDict{Symbol,Matrix{Float64}}()
    edges = nothing
    for (row, i) in enumerate(draws)
        hists = get_best_fit_histograms(samples.v[i], model, bin_edges)
        isnothing(groups) || (hists = group_histograms(hists, groups))
        for (name, h) in hists
            counts = get!(() -> zeros(n, length(h.weights)), predictions, name)
            counts[row, :] = h.weights
        end
        edges = first(values(hists)).edges
    end

    bands = OrderedDict{Symbol,NamedTuple{(:low, :high),NTuple{2,Histogram{Float64,1}}}}()
    for (name, counts) in predictions
        low, high = ([quantile(column, q) for column in eachcol(counts)] for q in quantiles)
        bands[name] = (low = Histogram(edges, low), high = Histogram(edges, high))
    end
    return bands
end

"""
    normalised_poisson_residual(mu, obs)

Residual of the observed counts `obs` with respect to the predicted counts `mu`,
in units of standard deviations. Takes two numbers and returns a `Float64`, or two
arrays of the same size and returns a `Vector{Float64}`.

For `mu >= 50` the Poisson distribution is close to a Gaussian and the residual is
`(obs - mu) / sqrt(mu)`. For smaller `mu` it is the p-value of `obs` (the probability
of at least `obs` counts above the mode, or of at most `obs` counts below it),
converted into the number of standard deviations of a Gaussian. The sign is positive
for an excess of data. A bin with `obs < 1` and a prediction below 1 has residual 0.
"""
function normalised_poisson_residual(mu::Real, obs::Real)
    x = obs < 1 ? 0 : obs
    if mu < 50
        mode = floor(Int, mu)
        x == 0 && mode == 0 && return 0.0

        if x < mode
            prob = Distributions.cdf(Distributions.Poisson(mu), x)
            sgn = -1
        else
            prob = 1 - Distributions.cdf(Distributions.Poisson(mu), x - 1)
            sgn = 1
        end
        return sgn * Distributions.quantile(Distributions.Normal(), 1 - prob)
    else
        return (x - mu) / sqrt(mu)
    end
end

function normalised_poisson_residual(mu::AbstractArray{<:Real}, obs::AbstractArray{<:Real})
    size(mu) == size(obs) || throw(DimensionMismatch("mu and obs differ in size"))
    return vec(normalised_poisson_residual.(mu, obs))
end
