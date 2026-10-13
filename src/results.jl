"""
    get_best_fit_histograms(parameters, model, bin_edges)

Histograms of the prediction of the model for the `parameters` (a named tuple, for
example the best fit). `model` and `bin_edges` come from
[`build_forward_model`](@ref).

Returns an ordered dictionary `name => Histogram` with the total prediction in `:total`,
followed by the contribution of each component under its name in the model, sorted by
name. The blinded region and
the range outside of the fit are empty bins (see [`counts_to_histogram`](@ref)).
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
