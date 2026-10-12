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
