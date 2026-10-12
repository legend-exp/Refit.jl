using IntervalSets
using Refit
using StatsBase
using Test

@testset "priors" begin
    edges = 0.0:100.0:1000.0
    pdf = Histogram(edges, fill(10.0, 10))
    # the data are lower in the bins 400-600 keV, where the peak is
    data_weights = fill(50.0, 10)
    data_weights[5:6] .= 20.0
    data = Histogram(edges, data_weights)

    @testset "get_max_norm" begin
        # 8 bins between 100 and 900 keV: 340 data counts, 80 simulated counts
        @test get_max_norm(100, 900, data, pdf) ≈ 340 / 80
        @test get_max_norm(100, 900, data, pdf; factor = 1.1, norm = 2) ≈ 1.1 * 340 / 80 / 2
        # a window edge that is not a bin edge: only the bins fully inside are used
        @test get_max_norm(420, 580, data, pdf) == Inf
        # blinded bins are not used
        @test get_max_norm(0, 1000, data, pdf; blind = (400, 600)) ≈ 400 / 80
        @test get_max_norm(400, 600, data, pdf; blind = (400, 600)) == Inf
        @test_throws ErrorException get_max_norm(
            0,
            1000,
            Histogram(0.0:200.0:1000.0, fill(1.0, 5)),
            pdf,
        )
    end

    @testset "get_priors" begin
        pdfs = Dict(:K40_hpge_bulk => pdf, :other_source => pdf)
        range_max = 1.1 * 340 / 80
        peak_max = 1.5 * 40 / 20

        # without regions only the fit range is used
        priors = get_priors(pdfs, data; e_low = 100, e_high = 900)
        @test priors isa Dict{Symbol,<:Interval}
        @test priors[:other_source] == 0.0 .. range_max
        @test priors[:K40_hpge_bulk] == 0.0 .. range_max

        # a peak: the window is 400 to 600 keV and gives the smaller limit
        ref = Dict(:K40 => 500)
        priors = get_priors(
            pdfs,
            data;
            ref_regions = ref,
            e_low = 100,
            e_high = 900,
            peak_width = 100,
        )
        @test priors[:K40_hpge_bulk] == 0.0 .. peak_max
        @test priors[:other_source] == 0.0 .. range_max

        # a tuple is used as the window, the longest matching key is used
        ref = Dict(:K4 => (0, 200), :K40 => (400, 600), :other => 500)
        priors = get_priors(
            pdfs,
            data;
            ref_regions = ref,
            e_low = 100,
            e_high = 900,
            peak_width = 100,
        )
        @test priors[:K40_hpge_bulk] == 0.0 .. peak_max
        @test priors[:other_source] == 0.0 .. peak_max

        # settings
        priors =
            get_priors(pdfs, data; e_low = 100, e_high = 900, range_factor = 1, norm = 2)
        @test priors[:other_source] == 0.0 .. (340 / 80 / 2)

        # the result can be changed by hand
        priors[:other_source] = 0.0 .. 3.0
        @test priors[:other_source] == 0.0 .. 3.0

        @test_throws ArgumentError get_priors(
            Dict(:empty => Histogram(edges, zeros(10))),
            data,
        )
        @test_throws ArgumentError get_priors(
            pdfs,
            data;
            ref_regions = Dict(:K40 => "peak"),
        )
    end
end
