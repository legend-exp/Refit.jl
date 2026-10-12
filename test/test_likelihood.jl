using Refit
using Test

@testset "Poisson Likelihood" begin
    obs = [10, 20, 30]
    pred = [12.0, 18.0, 29.0]
    @test poisson_likelihood(obs, pred) ≈ sum(@. obs * log(pred) - pred)
    @test_throws DimensionMismatch poisson_likelihood([1, 2, 3], [1.0, 2.0])
    @test isfinite(poisson_likelihood([0], [0.0]))
end
