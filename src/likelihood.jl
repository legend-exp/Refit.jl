"""
    poisson_likelihood(obs, pred)

Poisson log-likelihood of observed counts `obs` given predicted counts `pred`,
dropping the constant `-log(k!)` term. A tiny `eps()` is added to the prediction inside
the logarithm, so that a bin predicted empty gives a finite value. A negative prediction
has no Poisson probability and gives `-Inf`, so that a sampler or an optimiser rejects
these parameters.
"""
function poisson_likelihood(obs::AbstractVector, pred::AbstractVector)
    length(obs) == length(pred) || throw(DimensionMismatch("obs and pred differ in length"))
    # the number type of the prediction, so that derivatives (ForwardDiff) work
    ll = zero(promote_type(eltype(obs), eltype(pred)))
    for (k, λ) in zip(obs, pred)
        λ < 0 && return oftype(ll, -Inf)
        ll += k * log(λ + eps()) - λ
    end
    return ll
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
