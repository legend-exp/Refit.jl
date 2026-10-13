using LegendHDF5IO
using Refit
using StatsBase
using IntervalSets
using Test
using YAML

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

@testset "write best fit" begin
    edges = 0.0:100.0:1000.0
    pdfs =
        Dict(:a => Histogram(edges, fill(10.0, 10)), :b => Histogram(edges, fill(5.0, 10)))
    model, bin_edges = build_forward_model(
        pdfs;
        e_low = 100,
        e_high = 900,
        blind = (450.0, 550.0),
        norm = 1e-3,
    )
    hists = get_best_fit_histograms((a = 2.0, b = 4.0), model, bin_edges)

    # write to lh5 and read back
    file = tempname() * ".lh5"
    data = counts_to_histogram(fill(3.0, length(bin_edges.low)), bin_edges)
    write_best_fit(file, hists; data)
    @test read_hist("best_fit/total", file).weights == hists[:total].weights
    @test read_hist("best_fit/a", file).weights == hists[:a].weights
    @test read_hist("best_fit/b", file).edges == hists[:b].edges
    @test read_hist("data", file).weights == data.weights

    # without data nothing is written for it
    write_best_fit(file, hists)
    @test !lh5open(f -> haskey(f, "data"), file, "r")
end

@testset "activities" begin
    parameters = (a = 60.0, b = 160.0)
    intervals = (a = 40.0 .. 70.0, b = 150.0 .. 190.0)

    file = tempname() * ".yaml"
    write_activities(file, parameters, intervals)
    result = YAML.load_file(file)
    @test result["credible_level"] == 0.68
    @test Set(keys(result["activities"])) == Set(["a", "b"])
    @test result["activities"]["a"] ==
          Dict("best_fit" => 60.0, "low" => 40.0, "high" => 70.0)
    @test result["activities"]["b"] ==
          Dict("best_fit" => 160.0, "low" => 150.0, "high" => 190.0)

    write_activities(file, parameters, intervals; level = 0.95)
    @test YAML.load_file(file)["credible_level"] == 0.95
end
