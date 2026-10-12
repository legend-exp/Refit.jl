using DensityInterface
using Refit
using StatsBase
using Test

@testset "Poisson Likelihood" begin
    obs = [10, 20, 30]
    pred = [12.0, 18.0, 29.0]
    @test poisson_likelihood(obs, pred) ≈ sum(@. obs * log(pred) - pred)
    @test_throws DimensionMismatch poisson_likelihood([1, 2, 3], [1.0, 2.0])
    @test isfinite(poisson_likelihood([0], [0.0]))
end

@testset "get_counts" begin
    data = Histogram(0.0:1000.0:5000.0, [1.0, 2.0, 3.0, 4.0, 5.0])
    @test get_counts(data, (low = [1000.0, 2000.0], high = [2000.0, 3000.0])) == [2.0, 3.0]
    # model bins that merge two data bins
    @test get_counts(data, (low = [0.0, 2000.0], high = [2000.0, 4000.0])) == [3.0, 7.0]
    @test get_counts(data, (low = [0.0, 3000.0], high = [1000.0, 5000.0])) == [1.0, 9.0]
    @test get_counts(data, (low = [1000.0], high = [2000.0])) isa Vector{Float64}
    @test_throws ArgumentError get_counts(data, (low = [500.0], high = [2000.0]))
    @test_throws DimensionMismatch get_counts(data, (low = [0.0, 1000.0], high = [1000.0]))
end

@testset "build_likelihood" begin
    edges = 0.0:1000.0:5000.0
    pdfs = Dict(:a => Histogram(edges, [1.0, 2.0, 3.0, 4.0, 5.0]))
    model, bin_edges = build_forward_model(pdfs; e_low = 1000, e_high = 4000)
    data = Histogram(edges, [7.0, 2.0, 3.0, 5.0, 1.0])
    counts = get_counts(data, bin_edges)
    @test counts == [2.0, 3.0, 5.0]

    likelihood = build_likelihood(counts, model)
    params = (activities = (a = 1000.0,),)
    @test logdensityof(likelihood, params) ≈ poisson_likelihood(counts, model(params))
    @inferred logdensityof(likelihood, params)

    wrong_length = build_likelihood([1.0, 2.0], model)
    @test_throws DimensionMismatch logdensityof(wrong_length, params)
end
