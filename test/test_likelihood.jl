using Refit
using Test

@testset "Poisson Likelihood Tests" begin
    obs = [10, 20, 30]
    pred = [12.0, 18.0, 29.0]
    ll = poisson_likelihood(obs, pred)
    @test isfinite(ll)
end
