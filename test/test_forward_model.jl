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

    # hot-loop requirements: inferred types and a single allocation per call
    p = (a = 1.0, b = 0.5)
    @inferred fm(p)
    fm(p)
    @test (@allocated fm(p)) <= @allocated zeros(3)

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
        # more than 30 spectra: still one allocation per call, parameters in any order
        names = [Symbol("s$i") for i = 1:60]
        many = Dict(name => h1 for name in names)
        fm, _ = build_forward_model(many; e_low = 0, e_high = 5000)
        p = NamedTuple{Tuple(names)}(ntuple(i -> Float64(i), 60))
        shuffled = NamedTuple{Tuple(reverse(names))}(Tuple(reverse(values(p))))
        @test fm(p) ≈ sum(1:60) .* h1.weights
        @test fm(shuffled) ≈ fm(p)
        @test_throws FieldError fm(Base.structdiff(p, (s1 = 0,)))
        fm(p)
        @test (@allocated fm(p)) <= @allocated zeros(5)
    end
end
