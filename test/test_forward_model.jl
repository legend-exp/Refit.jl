using Refit
using StatsBase
using Test

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

    _, e_blind = build_forward_model(pdfs; e_low = 0, e_high = 5000, blind = (2000, 3000))
    @test e_blind.low == [0.0, 1000.0, 3000.0, 4000.0]
    @test e_blind.high == [1000.0, 2000.0, 4000.0, 5000.0]

    _, e_rebin = build_forward_model(pdfs; e_low = 0, e_high = 4000, binning = 2)
    @test e_rebin.low == [0.0, 2000.0]
    @test e_rebin.high == [2000.0, 4000.0]

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

        @test_throws ArgumentError build_forward_model(Dict(:a => (pdf = h1, nrm = 1)))
        @test_throws ArgumentError build_forward_model(Dict(:a => (activity = p -> p.x,)))
        @test_throws ArgumentError build_forward_model(Dict(:a => NamedTuple()))
    end
end
