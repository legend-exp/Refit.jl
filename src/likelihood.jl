
"""
    poisson_likelihood(obs::AbstractVector, pred::AbstractVector) -> Real

Optimised implementation of the Poisson likelihood function for observed and predicted counts.
- `obs`: Vector of observed counts (non-negative integers).
- `pred`: Vector of predicted counts (non-negative reals).
Returns the log-likelihood value.
"""
@inline function poisson_likelihood(obs::AbstractVector, pred::AbstractVector)
    ll = 0.0
    e0 = eps(ll)

    @fastmath @inbounds @simd for i in eachindex(obs)
        λ = pred[i] + e0      # avoid zero
        k = obs[i]
        ll += k*log(λ) - λ
    end
    return ll
end
