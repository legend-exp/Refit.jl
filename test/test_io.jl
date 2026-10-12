using Refit
using Test
using StatsBase

@testset "test_read_histograms" begin

    path = joinpath(@__DIR__, "test_files")
    file = path*"/test_histograms.lh5"
    @test read_hist("hist", file) isa Histogram
end
