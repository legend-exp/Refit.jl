"""
    poisson_likelihood(obs, pred)

Poisson log-likelihood of observed counts `obs` given predicted counts `pred`,
dropping the constant `-log(k!)` term. `pred` must be positive; a tiny `eps()` is
added to avoid `log(0)`.
"""
function poisson_likelihood(obs::AbstractVector, pred::AbstractVector)
    length(obs) == length(pred) || throw(DimensionMismatch("obs and pred differ in length"))
    return sum(k * log(λ + eps()) - λ for (k, λ) in zip(obs, pred))
end

"""
    build_likelihood(counts, model)

Build the Poisson log-likelihood of the measured `counts` as a function of the
parameters of `model`.

`counts` is a flat vector with one entry per bin, in the same order as the output
of `model(params)`. Get it from a histogram with [`get_counts`](@ref), so that the
bins are the ones the model predicts. `model` comes from
[`build_forward_model`](@ref).

Returns a log-density function: evaluate it with
`DensityInterface.logdensityof(likelihood, params)`.
"""
function build_likelihood(counts::AbstractVector{<:Real}, model::Function)
    likelihood(params::NamedTuple) = poisson_likelihood(counts, model(params))
    return logfuncdensity(likelihood)
end
