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

    # no normalisation: the raw histogram
    @test read_pdf(file; n_prim_field = nothing).weights == weights
    @test read_pdf(file; n_prim_field = nothing).edges ==
          read_hist("pdf/mul_surv", file).edges

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

@testset "read_pdf with several files" begin
    source = joinpath(@__DIR__, "test_files", "test_histograms.lh5")
    binning, weights = lh5open(source, "r") do f
        f["hist/binning/axis_0"], f["hist/weights"][:]
    end

    # two files with a different number of primaries
    files = [tempname() * ".lh5" for _ = 1:2]
    for (file, scale, n_prim) in zip(files, [1, 2], [1000, 500])
        lh5open(file, "w") do f
            f["pdf/mul_surv/binning/axis_0"] = binning
            f["pdf/mul_surv/weights"] = scale .* weights
            f["nr_sim_events"] = n_prim
        end
    end

    # the weights 3 and 1 are normalised to 0.75 and 0.25
    pdf = read_pdf(Dict(files[1] => 3.0, files[2] => 1.0))
    @test pdf isa Histogram
    @test pdf.weights ≈ (0.75 .* weights ./ 1000) .+ (0.25 .* 2 .* weights ./ 500)

    # integer weights, and the same result if all the weights are scaled
    @test read_pdf(Dict(files[1] => 3, files[2] => 1)).weights ≈ pdf.weights
    @test read_pdf(Dict(files[1] => 30.0, files[2] => 10.0)).weights ≈ pdf.weights

    # one file is the single-file result
    @test read_pdf(Dict(files[2] => 5.0)).weights ≈ read_pdf(files[2]).weights

    # no division by the number of primaries
    pdf = read_pdf(Dict(files[1] => 1.0, files[2] => 1.0); n_prim_field = nothing)
    @test pdf.weights ≈ 1.5 .* weights
end
