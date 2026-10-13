using LegendHDF5IO
using Refit
using StatsBase
using Test

@testset "test_read_histograms" begin

    path = joinpath(@__DIR__, "test_files")
    file = path*"/test_histograms.lh5"
    @test read_hist("hist", file) isa Histogram
end

@testset "read_pdf" begin
    # a file with two histograms and the number of simulated primaries
    source = joinpath(@__DIR__, "test_files", "test_histograms.lh5")
    binning, weights = lh5open(source, "r") do f
        f["hist/binning/axis_0"], f["hist/weights"][:]
    end

    file = tempname() * ".lh5"
    lh5open(file, "w") do f
        f["pdf/mul_surv/binning/axis_0"] = binning
        f["pdf/mul_surv/weights"] = weights
        f["pdf/mul_other/binning/axis_0"] = binning
        f["pdf/mul_other/weights"] = 2 .* weights
        f["nr_sim_events"] = 1000
        f["n_events"] = 500
    end

    @test Refit.read_scalar("nr_sim_events", file) == 1000

    pdf = read_pdf(file)
    @test pdf isa Histogram
    @test pdf.edges == read_hist("hist", source).edges
    @test pdf.weights ≈ weights ./ 1000

    # one group as a string or a vector
    @test read_pdf(file; groups = "pdf/mul_other").weights ≈ 2 .* weights ./ 1000
    @test read_pdf(file; groups = ["pdf/mul_other"]).weights ≈ 2 .* weights ./ 1000

    # groups are summed before the normalisation
    pdf = read_pdf(file; groups = ["pdf/mul_surv", "pdf/mul_other"])
    @test pdf.weights ≈ 3 .* weights ./ 1000

    # weights multiply each group before the sum
    pdf = read_pdf(file; groups = ["pdf/mul_surv", "pdf/mul_other"], weights = [0.7, 0.3])
    @test pdf.weights ≈ (0.7 .+ 0.3 .* 2) .* weights ./ 1000
    @test read_pdf(file; groups = "pdf/mul_other", weights = [0.5]).weights ≈
          0.5 .* 2 .* weights ./ 1000
    @test_throws DimensionMismatch read_pdf(
        file;
        groups = ["pdf/mul_surv", "pdf/mul_other"],
        weights = [1.0],
    )

    pdf = read_pdf(file; groups = "pdf/mul_other", n_prim_field = "n_events")
    @test pdf.weights ≈ 2 .* weights ./ 500

    # the file is not changed by scaling the result
    pdf.weights ./= 2
    @test read_hist("pdf/mul_other", file).weights == 2 .* weights
end
