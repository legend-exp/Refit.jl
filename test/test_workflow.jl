using DensityInterface
using Refit
using StatsBase
using Test

# The steps of the example in the README, on simulated spectra read from a file
@testset "fit workflow" begin
    file = joinpath(@__DIR__, "test_files", "test_histograms.lh5")
    h = read_hist("hist", file)

    mc_dict = Dict(
        :Ar39_liquid_argon => h,
        :two_vbb_hpge_bulk => Histogram(h.edges, reverse(h.weights)),
    )

    # measured spectrum: the model for known activities, in 10 keV bins
    true_params = (activities = (Ar39_liquid_argon = 2.0, two_vbb_hpge_bulk = 0.5),)
    data = Histogram(h.edges, 2.0 .* h.weights .+ 0.5 .* reverse(h.weights))

    model, bin_edges = build_forward_model(
        mc_dict;
        e_low = 100,
        e_high = 900,
        blind = (400, 500),
        binning = 2,
    )
    counts = get_counts(data, bin_edges)
    @test length(counts) == length(model(true_params))
    @test counts ≈ model(true_params)

    # back to a histogram: the fit range, with 0 counts in the blinded region
    fit = counts_to_histogram(model(true_params), bin_edges)
    @test first(fit.edges[1]) == 100 && last(fit.edges[1]) == 900
    @test sum(fit.weights) ≈ sum(counts)
    gap = findfirst(==(400), fit.edges[1])
    @test fit.edges[1][gap+1] == 500 && fit.weights[gap] == 0
    @test get_counts(fit, bin_edges) ≈ counts

    likelihood = build_likelihood(counts, model)
    logL(activities) = logdensityof(likelihood, (activities = activities,))

    # the likelihood is largest at the true activities
    at_truth = logL(true_params.activities)
    @test at_truth > logL((Ar39_liquid_argon = 2.2, two_vbb_hpge_bulk = 0.5))
    @test at_truth > logL((Ar39_liquid_argon = 2.0, two_vbb_hpge_bulk = 0.4))
end
