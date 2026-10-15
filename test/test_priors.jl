using DensityInterface
using Distributions
using IntervalSets
using Random
using Refit
using StatsBase
using Test
using ValueShapes

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
        # a window edge that is not a bin edge: only the bins fully inside are used
        @test get_max_norm(420, 580, data, pdf) == Inf
        # a region left out between the ranges of the binning, such as a blinded one,
        # is not used
        blinded = [0:100:400, 600:100:1000]
        @test get_max_norm(0, 1000, data, pdf; binning = blinded) ≈ 400 / 80
        @test get_max_norm(400, 600, data, pdf; binning = blinded) == Inf
        # a binning without ranges uses all bins
        @test get_max_norm(100, 900, data, pdf; binning = 2) ≈ 340 / 80
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
        @test all(p -> p isa Interval, values(priors))
        @test priors[:other_source] == 0.0 .. range_max
        @test priors[:K40_hpge_bulk] == 0.0 .. range_max

        # a reference range of 400 to 600 keV gives the smaller limit
        ref = Dict(:K40 => (400, 600))
        priors = get_priors(pdfs, data; ref_ranges = ref, e_low = 100, e_high = 900)
        @test priors[:K40_hpge_bulk] == 0.0 .. peak_max
        @test priors[:other_source] == 0.0 .. range_max

        # the longest matching key is used
        ref = Dict(:K4 => (0, 200), :K40 => (400, 600), :other => (400.0, 600.0))
        priors = get_priors(pdfs, data; ref_ranges = ref, e_low = 100, e_high = 900)
        @test priors[:K40_hpge_bulk] == 0.0 .. peak_max
        @test priors[:other_source] == 0.0 .. peak_max

        # safety factors
        priors = get_priors(
            pdfs,
            data;
            ref_ranges = Dict(:K40 => (400, 600)),
            e_low = 100,
            e_high = 900,
            range_factor = 2,
            peak_factor = 3,
        )
        @test priors[:K40_hpge_bulk] == 0.0 .. (3 * 40 / 20)
        @test priors[:other_source] == 0.0 .. (2 * 340 / 80)

        # norm, and a blinded region left out between two ranges
        priors = get_priors(pdfs, data; e_low = 100, e_high = 900, norm = 2)
        @test priors[:other_source] == 0.0 .. (range_max / 2)
        priors = get_priors(
            pdfs,
            data;
            ref_ranges = Dict(:K40 => (400, 600)),
            e_low = 100,
            e_high = 900,
            binning = [100:100:400, 600:100:900],
        )
        @test priors[:K40_hpge_bulk] == 0.0 .. (1.1 * 300 / 60)

        # the result can be changed by hand
        priors[:other_source] = 0.0 .. 3.0
        @test priors[:other_source] == 0.0 .. 3.0

        @test_throws ArgumentError get_priors(
            Dict(:empty => Histogram(edges, zeros(10))),
            data,
        )
        @test_throws TypeError get_priors(pdfs, data; ref_ranges = Dict(:K40 => 1461))

        # a generalised histogram: the widest range of its grid points, here the one
        # with the fewest counts
        gh = GeneralisedHistogram([pdf, Histogram(edges, fill(5.0, 10))], :x => 0.0:1.0:1.0)
        priors = get_priors(Dict(:g => gh), data; e_low = 100, e_high = 900)
        @test priors[:g] == 0.0 .. (1.1 * 340 / 40)
    end
end

# The ranges become a prior for the fit. BAT samples a prior made of ValueShapes
# distributions, which needs the parameters to be a flat named tuple.
@testset "prior distribution" begin
    edges = 0.0:100.0:1000.0
    pdfs =
        Dict(:a => Histogram(edges, fill(10.0, 10)), :b => Histogram(edges, fill(5.0, 10)))
    data = Histogram(edges, fill(50.0, 10))

    model, bin_edges = build_forward_model(pdfs; e_low = 100, e_high = 900, norm = 1e-3)
    counts = get_counts(data, bin_edges)
    likelihood = build_likelihood(counts, model)
    ranges = get_priors(pdfs, data; e_low = 100, e_high = 900, norm = 1e-3)

    # components with an activity or analytic function are skipped, the others are as
    # without them
    mixed = merge(
        pdfs,
        Dict(
            :c => (pdf = pdfs[:a], activity = p -> p.s),
            :d => (analytic = (E, p) -> p.r,),
        ),
    )
    # in any order, the dictionary has none
    mixed_ranges = @test_logs (:info, r"c has an activity function") (
        :info,
        r"d has an analytic function",
    ) match_mode = :any get_priors(mixed, data; e_low = 100, e_high = 900, norm = 1e-3)
    @test mixed_ranges == ranges

    # a uniform distribution on each range
    prior = build_prior_dist(ranges)
    @test prior isa NamedTupleDist

    # a distribution is used as it is
    ranges[:K40_hpge_bulk] = truncated(Normal(50.0, 10.0), 0.0, Inf)
    @test build_prior_dist(ranges) isa NamedTupleDist

    Random.seed!(1)  # fixed seed, so that the random draw is the same in every run
    draw = rand(prior)
    @test draw isa NamedTuple
    @test Set(keys(draw)) == Set(keys(pdfs))
    @test all(ranges[k].left <= draw[k] <= ranges[k].right for k in keys(draw))

    # a draw goes directly into the model and the likelihood, in any order
    @test length(model(draw)) == length(counts)
    @test isfinite(logdensityof(likelihood, draw))
    reversed = NamedTuple{reverse(keys(draw))}(draw)
    @test logdensityof(likelihood, reversed) == logdensityof(likelihood, draw)

    # nested named tuples cannot be sampled, so the activities are not nested
    @test_throws MethodError rand(NamedTupleDist(activities = prior))
end
