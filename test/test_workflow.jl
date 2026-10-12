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

    # back to a histogram: equal to the data in the fit range, 0 outside of it
    fit = counts_to_histogram(model(true_params), bin_edges, rebin_hist(data, 2))
    in_fit = Refit.get_valid_bins(fit.edges[1], (100, 900), (400, 500))
    @test fit.weights[in_fit] ≈ rebin_hist(data, 2).weights[in_fit]
    @test all(fit.weights[.!in_fit] .== 0)

    likelihood = build_likelihood(counts, model)
    logL(activities) = logdensityof(likelihood, (activities = activities,))

    # the likelihood is largest at the true activities
    at_truth = logL(true_params.activities)
    @test at_truth > logL((Ar39_liquid_argon = 2.2, two_vbb_hpge_bulk = 0.5))
    @test at_truth > logL((Ar39_liquid_argon = 2.0, two_vbb_hpge_bulk = 0.4))
end
