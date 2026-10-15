using Refit
using Test
using StatsBase
using Random

# a histogram with the given counts, in 1 keV bins from 0 keV
hist_with(weights) = Histogram((0.0:1.0:length(weights),), Float64.(weights))

@testset "GeneralisedHistogram" begin

    @testset "one parameter" begin
        hists = [hist_with([1, 2, 3]), hist_with([3, 2, 1])]
        gh = GeneralisedHistogram(hists; fccd = 1.0:0.5:1.5)
        @test gh.grid == (fccd = 1.0:0.5:1.5,)

        # the grid points themselves, and halfway
        @test get_weights(gh, (fccd = 1.0,)) == [1, 2, 3]
        @test get_weights(gh, (fccd = 1.5,)) == [3, 2, 1]
        @test get_weights(gh, (fccd = 1.25,)) ≈ [2, 2, 2]
        # other parameters are ignored
        @test get_weights(gh, (A = 10.0, fccd = 1.25)) ≈ [2, 2, 2]

        h = get_histogram(gh, (fccd = 1.25,))
        @test h isa Histogram
        @test h.edges[1] == 0:1:3
        @test h.weights ≈ [2, 2, 2]

        # outside the grid
        @test_throws DomainError get_weights(gh, (fccd = 0.9,))
        @test_throws DomainError get_weights(gh, (fccd = 1.6,))
    end

    # multilinear interpolation reproduces a function that is linear in each parameter
    # exactly, also between grid points
    @testset "$N parameters" for N = 1:6
        names = (:a, :b, :c, :d, :e, :f)[1:N]
        grid = NamedTuple{names}(ntuple(d -> range(0.0, step = 0.5 * d, length = 3), N))
        f(x) = 1 + sum(d * x[d] for d = 1:N) + (N > 1 ? x[1] * x[2] : 0.0)
        hists = map(CartesianIndices(map(length, Tuple(grid)))) do idx
            x = ntuple(d -> grid[d][idx[d]], N)
            hist_with([f(x), 2 * f(x)])
        end
        gh = GeneralisedHistogram(hists; grid...)
        # the same from a function of the parameter values
        spectrum(p) = (x = Tuple(p); hist_with([f(x), 2 * f(x)]))
        @test GeneralisedHistogram(spectrum; grid...).counts == gh.counts

        rng = MersenneTwister(N)
        for _ = 1:20
            x = ntuple(d -> last(grid[d]) * rand(rng), N)
            @test get_weights(gh, NamedTuple{names}(x)) ≈ [f(x), 2 * f(x)]
        end
        # the order of the parameters does not matter
        x = ntuple(d -> last(grid[d]) / 3, N)
        @test get_weights(gh, NamedTuple{reverse(names)}(reverse(x))) ≈ [f(x), 2 * f(x)]
        # all types are known when Julia compiles the function, which makes it fast
        @test (@inferred get_weights(gh, NamedTuple{names}(x))) isa Vector{Float64}
    end

    @testset "rebinning" begin
        rng = MersenneTwister(1)
        hists = [hist_with(rand(rng, 12)) for _ = 1:3, _ = 1:4]
        gh = GeneralisedHistogram(hists; z = -1.0:1.0:1.0, phi = 0.0:2.0:6.0)
        params = (z = 0.3, phi = 4.5)

        # rebinning before or after the interpolation gives the same counts
        for binning in (3, [0.0, 2.0, 7.0, 12.0])
            r = rebin_hist(gh, binning)
            @test r.grid == gh.grid
            @test r.edges == rebin_hist(get_histogram(gh, params), binning).edges[1]
            @test get_weights(r, params) ≈
                  rebin_hist(get_histogram(gh, params), binning).weights
        end
    end

    @testset "add_histograms" begin
        gh1 = GeneralisedHistogram([hist_with([1, 2]), hist_with([3, 4])]; x = 0.0:1.0:1.0)
        gh2 = GeneralisedHistogram([hist_with([5, 6]), hist_with([7, 8])]; x = 0.0:1.0:1.0)
        total = add_histograms([gh1, gh2])
        @test total.counts == [[6, 8], [10, 12]]
        @test get_weights(total, (x = 0.5,)) ≈ [8, 10]

        gh3 = GeneralisedHistogram([hist_with([5, 6]), hist_with([7, 8])]; x = 0.0:2.0:2.0)
        @test_throws ErrorException add_histograms([gh1, gh3])
    end

    @testset "invalid input" begin
        two = [hist_with([1, 2]), hist_with([3, 4])]
        # different binnings
        bad = [hist_with([1, 2]), hist_with([1, 2, 3])]
        @test_throws ArgumentError GeneralisedHistogram(bad; x = 0.0:1.0:1.0)
        # the grid does not match the array of histograms
        @test_throws DimensionMismatch GeneralisedHistogram(two; x = 0.0:1.0:2.0)
        @test_throws ArgumentError GeneralisedHistogram(two; x = 0.0:1.0:1.0, y = 0:1)
        # the grid is not an evenly spaced range of at least two values
        @test_throws ArgumentError GeneralisedHistogram(two; x = [0.0, 1.0])
        @test_throws ArgumentError GeneralisedHistogram(two[1:1]; x = 0.0:1.0:0.0)
        # no parameter
        @test_throws ArgumentError GeneralisedHistogram(two)
        # an array with its dimensions in the wrong order, when the grids differ in length
        hists = [hist_with([x, y]) for y = 0:2, x = 0:1]
        @test_throws DimensionMismatch GeneralisedHistogram(hists; x = 0:1, y = 0:2)
    end
end
