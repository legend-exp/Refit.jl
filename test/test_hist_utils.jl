using Refit
using Test
using StatsBase

@testset "Histogram Utilities" begin
    h = fit(Histogram, randn(1000), -5:0.5:5)

    @testset "rebin_hist" begin
        r = rebin_hist(h, 2)
        @test length(r.weights) == 10
        @test sum(r.weights) == sum(h.weights)

        # not divisible, with range edges
        r = rebin_hist(h, 3)
        @test length(r.weights) == 7
        @test sum(r.weights) == sum(h.weights)
        @test r.edges[1][end] == 5.0
        @test_throws ArgumentError rebin_hist(h, 0)
    end

    @testset "add_histograms" begin
        h1 = fit(Histogram, randn(1000), -5:0.5:5)
        h2 = fit(Histogram, randn(1000), -5:0.5:5)
        @test add_histograms([h1, h2]).weights == h1.weights .+ h2.weights
        @test eltype(add_histograms([h1, h2]).weights) == Int
        @test_throws ErrorException add_histograms([h1, fit(Histogram, randn(10), -5:1:5)])
    end

    @testset "integrate_hist" begin
        hc = Histogram(0:0.5:2, [10, 20, 30, 40])
        @test integrate_hist(hc, 0, 1) == 30
        @test integrate_hist(hc, 0.5, 2) == 90
        @test integrate_hist(hc, 0, 2) == sum(hc.weights)
        @test_throws ArgumentError integrate_hist(hc, 0.25, 1)
        @test_throws ArgumentError integrate_hist(hc, 0, 10)
    end

    @testset "get_valid_bins" begin
        edges = [0.0, 1, 2, 3, 4]
        @test Refit.get_valid_bins(edges, (1, 4)) == [false, true, true, true]
        @test Refit.get_valid_bins(edges, (0, 4), (1, 3)) == [true, false, false, true]
    end
end
