using Refit
using StatsBase
using Test

# global variables read by analytic functions in the tests below: an ordinary one makes
# the model slow, a constant does not
slow_rate = 2.0
const const_rate = 2.0

@testset "build_forward_model" begin
    edges = 0.0:1000.0:5000.0
    h1 = Histogram(edges, [1.0, 2.0, 3.0, 4.0, 5.0])
    h2 = Histogram(edges, [10.0, 20.0, 30.0, 40.0, 50.0])
    pdfs = Dict(:a => h1, :b => h2)

    fm, bin_edges = build_forward_model(pdfs; e_low = 1000, e_high = 4000, norm = 2e-3)
    @test bin_edges.low == [1000.0, 2000.0, 3000.0]
    @test bin_edges.high == [2000.0, 3000.0, 4000.0]
    @test fm((a = 1.0, b = 0.0)) ≈ [2.0, 3.0, 4.0] .* 2e-3
    @test fm((a = 1.0, b = 0.5)) ≈ ([2.0, 3.0, 4.0] .+ 0.5 .* [20.0, 30.0, 40.0]) .* 2e-3
    @test_throws FieldError fm((a = 1.0,))
    @test fm((b = 0.5, a = 1.0)) ≈ fm((a = 1.0, b = 0.5))
    @inferred fm((b = 0.5, a = 1.0))

    # hot-loop requirements: inferred types and two arrays per call, the activities and
    # the prediction
    p = (a = 1.0, b = 0.5)
    @inferred fm(p)
    fm(p)
    @test (@allocated fm(p)) <= (@allocated zeros(2)) + (@allocated zeros(3))

    # a blinded region, left out between two ranges
    _, e_blind = build_forward_model(
        pdfs;
        e_low = 0,
        e_high = 5000,
        binning = [[0, 1000, 2000], [3000, 4000, 5000]],
    )
    @test e_blind.low == [0.0, 1000.0, 3000.0, 4000.0]
    @test e_blind.high == [1000.0, 2000.0, 4000.0, 5000.0]

    _, e_rebin = build_forward_model(pdfs; e_low = 0, e_high = 4000, binning = 2)
    @test e_rebin.low == [0.0, 2000.0]
    @test e_rebin.high == [2000.0, 4000.0]

    # bins of any width, given by their edges
    fm_var, e_var = build_forward_model(
        pdfs;
        e_low = 1000,
        e_high = 5000,
        binning = [0.0, 1000.0, 3000.0, 5000.0],
    )
    @test e_var.low == [1000.0, 3000.0]
    @test e_var.high == [3000.0, 5000.0]
    @test fm_var((a = 1.0, b = 0.0)) ≈ [2.0 + 3.0, 4.0 + 5.0]
    # an analytic component gets the width of each bin
    fm_flat, e_flat = build_forward_model(
        Dict(:c => (analytic = (E, p) -> p.rate,));
        e_low = 0,
        e_high = 5000,
        edges = 0.0:1000.0:5000.0,
        binning = [0.0, 2000.0, 5000.0],
    )
    @test e_flat.high == [2000.0, 5000.0]
    @test fm_flat((rate = 2.0,)) ≈ [4000.0, 6000.0]

    # several ranges, each with its own edges: the region between them is skipped
    ranges = [[0.0, 1000.0], [2000.0, 3000.0, 5000.0]]
    fm_r, e_r = build_forward_model(pdfs; e_low = 0, e_high = 5000, binning = ranges)
    @test e_r.low == [0.0, 2000.0, 3000.0]
    @test e_r.high == [1000.0, 3000.0, 5000.0]
    @test fm_r((a = 1.0, b = 0.0)) ≈ [1.0, 3.0, 4.0 + 5.0]
    # as read from a file, in any order, and with e_low still applied
    _, e_any = build_forward_model(
        pdfs;
        e_low = 500,
        e_high = 5000,
        binning = Any[Any[2000, 3000, 5000], Any[0, 1000]],
    )
    @test e_any.low == [2000.0, 3000.0]
    # ranges that touch leave no gap
    _, e_touch = build_forward_model(
        pdfs;
        e_low = 0,
        e_high = 5000,
        binning = [[0, 1000], [1000, 3000]],
    )
    @test e_touch.high == [1000.0, 3000.0]
    @test_throws ArgumentError build_forward_model(
        pdfs;
        binning = [[0, 3000], [2000, 5000]],
    )
    @test_throws ArgumentError build_forward_model(pdfs; binning = [[0, 1000], [2000]])

    # a spec with Julia ranges, or the same with strings as read from a YAML file
    h10 = Histogram(0.0:10.0:5000.0, ones(500))
    spec = [100:20:3000, 3200:100:4200, [4200, 4250, 4270]]
    _, e_spec =
        build_forward_model(Dict(:a => h10); e_low = 0, e_high = 5000, binning = spec)
    @test e_spec.low == [100.0:20.0:2980.0; 3200.0:100.0:4100.0; 4200.0; 4250.0]
    @test e_spec.high == [120.0:20.0:3000.0; 3300.0:100.0:4200.0; 4250.0; 4270.0]
    _, e_yaml = build_forward_model(
        Dict(:a => h10);
        e_low = 0,
        e_high = 5000,
        binning = Any["100:20:3000", "3200:100:4200", [4200, 4250, 4270]],
    )
    @test e_yaml == e_spec
    # a string ends at its stop, with a narrower last bin, so ranges that touch leave no
    # gap. A Julia range stops at its last step, so here 1900 to 2010 keV is skipped
    _, e_stop = build_forward_model(
        Dict(:a => h10);
        e_low = 0,
        e_high = 5000,
        binning = ["500:200:2010", "2010:10:2040"],
    )
    @test (1900.0, 2010.0) in collect(zip(e_stop.low, e_stop.high))
    @test e_stop.high[1:(end-1)] == e_stop.low[2:end]
    _, e_julia = build_forward_model(
        Dict(:a => h10);
        e_low = 0,
        e_high = 5000,
        binning = [500:200:2010, 2010:10:2040],
    )
    @test !((1900.0, 2010.0) in collect(zip(e_julia.low, e_julia.high)))
    # a decimal step that reaches its stop gives no extra bin
    _, e_dec = build_forward_model(
        Dict(:a => Histogram(0.0:0.1:1.0, ones(10)));
        e_low = 0,
        e_high = 1,
        binning = ["0:0.1:1"],
    )
    @test length(e_dec.low) == 10
    @test_throws ArgumentError build_forward_model(Dict(:a => h10); binning = ["100:3000"])
    @test_throws ArgumentError build_forward_model(
        Dict(:a => h10);
        binning = [100:20:200, :x],
    )

    @test_throws MethodError build_forward_model(Dict("a" => h1))
    @test_throws MethodError build_forward_model(Dict(:a => [1.0, 2.0]))
    @test_throws ArgumentError build_forward_model(Dict{Symbol,typeof(h1)}())

    @testset "many spectra" begin
        # more than 30 spectra: still two arrays per call, parameters in any order
        names = [Symbol("s$i") for i = 1:60]
        many = Dict(name => h1 for name in names)
        fm, _ = build_forward_model(many; e_low = 0, e_high = 5000)
        p = NamedTuple{Tuple(names)}(ntuple(i -> Float64(i), 60))
        shuffled = NamedTuple{Tuple(reverse(names))}(Tuple(reverse(values(p))))
        @test fm(p) ≈ sum(1:60) .* h1.weights
        @test fm(shuffled) ≈ fm(p)
        @test_throws FieldError fm(Base.structdiff(p, (s1 = 0,)))
        fm(p)
        @test (@allocated fm(p)) <= (@allocated zeros(60)) + (@allocated zeros(5))
    end

    @testset "activity functions" begin
        # :b is tied to the activity of :a, :c has its activity from a half-life
        comps = Dict(
            :a => h1,
            :b => (pdf = h2, activity = p -> 0.5 * p.a),
            :c => (pdf = h1, activity = p -> log(2) / p.t),
        )
        fm, _ = build_forward_model(comps; e_low = 1000, e_high = 4000, norm = 2e-3)
        p = (a = 1.0, t = 2.0)
        # the global norm applies to every component
        expected_a = [2.0, 3.0, 4.0] .* 2e-3
        expected_b = 0.5 .* [20.0, 30.0, 40.0] .* 2e-3
        expected_c = log(2) / 2.0 .* [2.0, 3.0, 4.0] .* 2e-3
        @test fm(p) ≈ expected_a .+ expected_b .+ expected_c

        parts = fm(p, Val(:by_component))
        @test keys(parts) == (:a, :b, :c)
        @test parts.a ≈ expected_a
        @test parts.b ≈ expected_b
        @test parts.c ≈ expected_c
        @test sum(values(parts)) ≈ fm(p)

        # a parameter of a function missing from `params` is an error
        @test_throws FieldError fm((a = 1.0,))

        # two arrays, plus less than 64 bytes for each of the two functions, to pass the
        # parameters to it
        @inferred fm(p)
        fm(p)
        @test (@allocated fm(p)) <= 2 * (@allocated zeros(3)) + 2 * 64

        # a pdf without an activity function is the same as giving the histogram
        fm_pdf, _ = build_forward_model(
            Dict(:a => (pdf = h1,));
            e_low = 1000,
            e_high = 4000,
            norm = 2e-3,
        )
        @test fm_pdf((a = 1.0,)) ≈ expected_a

        # histograms with integer counts work as well
        h_int = Histogram(edges, [1, 2, 3, 4, 5])
        fm_int, _ = build_forward_model(
            Dict(:a => h_int, :b => (pdf = h_int, activity = p -> 2 * p.a));
            e_low = 1000,
            e_high = 4000,
            norm = 2e-3,
        )
        @test fm_int((a = 1.0,)) ≈ 3 .* expected_a

        @test_throws ArgumentError build_forward_model(Dict(:a => (pdf = h1, nrm = 1)))
        @test_throws ArgumentError build_forward_model(Dict(:a => (activity = p -> p.x,)))
        @test_throws ArgumentError build_forward_model(Dict(:a => NamedTuple()))
    end
    @testset "analytic functions" begin
        # a histogram, a histogram times a function of the energy, an analytic component
        comps = Dict(
            :a => h1,
            :b => (pdf = h2, analytic = (E, p) -> p.eff * E / 1000),
            :c => (analytic = (E, p) -> p.rate * exp(-E / 1000),),
        )
        fm, _ = build_forward_model(comps; e_low = 1000, e_high = 4000, norm = 2e-3)
        centre = [1500.0, 2500.0, 3500.0]
        p = (a = 1.0, b = 2.0, eff = 0.5, rate = 3.0)
        expected_a = [2.0, 3.0, 4.0] .* 2e-3
        # :b has the activity p.b, the global norm and the function of the energy
        expected_b = 2.0 .* 2e-3 .* (0.5 .* centre ./ 1000) .* [20.0, 30.0, 40.0]
        # the global norm does not apply to :c, its function is in counts per keV
        expected_c = 3.0 .* exp.(-centre ./ 1000) .* 1000.0
        @test fm(p) ≈ expected_a .+ expected_b .+ expected_c

        parts = fm(p, Val(:by_component))
        @test keys(parts) == (:a, :b, :c)
        @test parts.a ≈ expected_a
        @test parts.b ≈ expected_b
        @test parts.c ≈ expected_c

        # a function of the energy together with an activity function
        fm_both, _ = build_forward_model(
            Dict(:b =>
                    (pdf = h2, activity = p -> 2 * p.x, analytic = (E, p) -> E / 1000));
            e_low = 1000,
            e_high = 4000,
            norm = 2e-3,
        )
        @test fm_both((x = 1.0,)) ≈ 2.0 .* 2e-3 .* centre ./ 1000 .* [20.0, 30.0, 40.0]

        # several points per bin. For g = E^2 the average over n points of a bin of width
        # w is below the true average by w^2 / (12 n^2), and the integral by w^3 / (12 n^2)
        square = Dict(:s => (analytic = (E, p) -> E^2,))
        for n in (1, 10)
            fm_n, _ = build_forward_model(
                square;
                e_low = 0,
                e_high = 3,
                edges = 0.0:1.0:3.0,
                points_per_bin = n,
            )
            @test fm_n((;)) ≈ [1 / 3, 7 / 3, 19 / 3] .- 1 / (12 * n^2)
        end
        # the same average multiplies a simulated spectrum
        fm_sq, _ = build_forward_model(
            Dict(:a => (pdf = h1, analytic = (E, p) -> E^2));
            e_low = 1000,
            e_high = 4000,
            norm = 2e-3,
            points_per_bin = 10,
        )
        true_mean = [7e9, 19e9, 37e9] ./ 3000  # of E^2 over each bin of 1000 keV
        @test fm_sq((a = 1.0,)) ≈ 2e-3 .* [2.0, 3.0, 4.0] .* (true_mean .- 1e6 / 1200)
        @test_throws ArgumentError build_forward_model(
            square;
            edges = 0.0:1.0:3.0,
            points_per_bin = 0,
        )

        # two arrays, plus a few bytes for each of the two functions, to pass them the
        # parameters. Not more in every bin
        @inferred fm(p)
        fm(p)
        @test (@allocated fm(p)) <= 2 * (@allocated zeros(3)) + 2 * 96

        # analytic components only: the binning comes from `edges`, rebinned as well
        only_g = Dict(:c => (analytic = (E, p) -> p.rate,))
        fm_g, e_g = build_forward_model(
            only_g;
            e_low = 1000,
            e_high = 4000,
            edges = 0.0:1000.0:5000.0,
        )
        @test e_g.low == [1000.0, 2000.0, 3000.0]
        @test fm_g((rate = 2.0,)) ≈ fill(2000.0, 3)
        _, e_g2 = build_forward_model(
            only_g;
            e_low = 0,
            e_high = 4000,
            edges = 0.0:1000.0:5000.0,
            binning = 2,
        )
        @test e_g2.low == [0.0, 2000.0]
        @test_throws ArgumentError build_forward_model(only_g)
        # `edges` must agree with the histograms
        @test_throws ErrorException build_forward_model(
            Dict(:a => h1);
            edges = 0.0:500.0:5000.0,
        )
        # a warning for an analytic function that reads an ordinary global variable
        @test_logs (:warn, r"slow_rate") build_forward_model(
            Dict(:c => (analytic = (E, p) -> slow_rate,));
            edges = 0.0:1000.0:5000.0,
        )
        @test_logs build_forward_model(
            Dict(:c => (analytic = (E, p) -> const_rate,));
            edges = 0.0:1000.0:5000.0,
        )
        # an analytic component with no pdf has no activity
        @test_throws ArgumentError build_forward_model(
            Dict(:c => (analytic = (E, p) -> p.rate, activity = p -> p.x));
            edges = 0.0:1000.0:5000.0,
        )
    end
end
