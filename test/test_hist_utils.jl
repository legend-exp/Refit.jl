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

        # to given edges, of any width, which must be edges of `h`
        edges = [-5.0, -1.0, 0.0, 0.5, 5.0]
        r = rebin_hist(h, edges)
        @test r.edges[1] == edges
        @test r.weights == [integrate_hist(h, edges[i], edges[i+1]) for i = 1:4]
        @test sum(r.weights) == sum(h.weights)
        # only the range of the edges
        @test rebin_hist(h, [-1.0, 1.0]).weights == [integrate_hist(h, -1.0, 1.0)]
        @test_throws ArgumentError rebin_hist(h, [-5.0, 0.25, 5.0])
        @test_throws ArgumentError rebin_hist(h, [0.0, -1.0])
        @test_throws ArgumentError rebin_hist(h, [0.0])
    end

    @testset "select_bins" begin
        h1 = Histogram(0.0:1.0:10.0, collect(1.0:10.0))
        # by default all bins
        counts, bin_edges = select_bins(h1)
        @test counts == 1:10
        @test bin_edges.low == 0:9 && bin_edges.high == 1:10
        # rebinned, then only the bins fully inside the range
        counts, bin_edges = select_bins(h1; binning = 2, e_low = 1, e_high = 9)
        @test counts == [3 + 4, 5 + 6, 7 + 8]
        @test bin_edges.low == [2, 4, 6] && bin_edges.high == [4, 6, 8]
        # a region left out between two ranges
        counts, bin_edges = select_bins(h1; binning = [0:2:4, 6:2:10])
        @test counts == [1 + 2, 3 + 4, 7 + 8, 9 + 10]
        @test bin_edges.low == [0, 2, 6, 8]
        # the same bins and counts as a forward model
        model, model_edges =
            build_forward_model(Dict(:a => h1); e_low = 1, e_high = 9, binning = 2)
        counts, bin_edges = select_bins(h1; binning = 2, e_low = 1, e_high = 9)
        @test bin_edges == model_edges
        @test counts == get_counts(h1, model_edges) == model((a = 1.0,))
    end

    @testset "uniform_binning_with_peaks" begin
        h1 = Histogram(0.0:1.0:100.0, ones(100))
        # no peaks: uniform bins, the last one narrower
        @test uniform_binning_with_peaks(h1, 30) == [0.0, 30.0, 60.0, 90.0, 100.0]
        # a peak from 44 to 52 is one bin. The edges 40, 50 and 60, less than one bin
        # width from it, go, so the bins next to it are 14 and 18 wide
        peak = [(44, 52)]
        @test uniform_binning_with_peaks(h1, 10; peaks = peak) ==
              [0.0, 10, 20, 30, 44, 52, 70, 80, 90, 100]
        # the edges are moved to the nearest edges of the histogram
        h2 = Histogram(0.0:2.0:100.0, ones(50))
        @test uniform_binning_with_peaks(h2, 10; peaks = [(44.6, 51.3)]) ==
              uniform_binning_with_peaks(h1, 10; peaks = peak)
        @test uniform_binning_with_peaks(h1, 10; low = 20, high = 50) == [20.0, 30, 40, 50]
        # usable as the binning of a model
        _, bin_edges = build_forward_model(
            Dict(:a => h1);
            e_low = 0,
            e_high = 100,
            binning = uniform_binning_with_peaks(h1, 10; peaks = peak),
        )
        @test (44.0, 52.0) in collect(zip(bin_edges.low, bin_edges.high))
        # a blinded region is a gap between two ranges, its neighbours are as for a peak
        blind = uniform_binning_with_peaks(h1, 10; peaks = peak, blinded = (74, 82))
        @test blind == [[0.0, 10, 20, 30, 44, 52, 74], [82.0, 100]]
        _, bin_edges =
            build_forward_model(Dict(:a => h1); e_low = 0, e_high = 100, binning = blind)
        @test all(@. bin_edges.high <= 74 || bin_edges.low >= 82)
        @test (82.0, 100.0) in collect(zip(bin_edges.low, bin_edges.high))
        # a blinded region at the end of the range leaves one range
        @test uniform_binning_with_peaks(h1, 10; blinded = (90, 100)) ==
              [[0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90]]
        @test_throws ArgumentError uniform_binning_with_peaks(
            h1,
            10;
            peaks = peak,
            blinded = (50, 60),
        )
        @test_throws ArgumentError uniform_binning_with_peaks(h1, 0)
        @test_throws ArgumentError uniform_binning_with_peaks(h1, 10; peaks = [(90, 120)])
        @test_throws ArgumentError uniform_binning_with_peaks(
            h1,
            10;
            peaks = [(40, 50), (45, 55)],
        )
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
        # only the bins fully inside the range
        @test Refit.get_valid_bins(edges, (0.5, 3)) == [false, true, true, false]
    end
end
