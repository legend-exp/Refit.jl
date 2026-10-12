using Refit
using Test
using StatsBase

@testset "Histogram Utilities" begin

    # Test rebin_hist
    h = fit(Histogram, randn(1000), -5:0.5:5)
    rebinned_h = Refit.rebin_hist(h, 2)
    @test length(rebinned_h.weights) == ceil(Int, length(h.weights) / 2)

    # Test histogram addition
    h1 = fit(Histogram, randn(1000), -5:0.5:5)
    h2 = fit(Histogram, randn(1000), -5:0.5:5)
    h_sum = Refit.add_histograms([h1, h2])
    @test h_sum.weights == h1.weights .+ h2.weights

    # Test integration
    a, b = -1.0, 1.0
    integral = Refit.integrate_hist(h, a, b)
    @test integral >= 0  # Integral should be non-negative
end
