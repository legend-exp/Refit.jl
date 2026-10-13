using LegendHDF5IO
using Refit
using StatsBase
using Test
using YAML

@testset "best fit histograms" begin
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
    parameters = (a = 2.0, b = 4.0)

    hists = get_best_fit_histograms(parameters, model, bin_edges)
    @test Set(keys(hists)) == Set([:total, :a, :b])
    @test all(h -> h isa Histogram, values(hists))

    # the sources add up to the total
    @test hists[:total].weights ≈ hists[:a].weights .+ hists[:b].weights
    @test hists[:a].weights ≈ model((a = 2.0, b = 0.0)) |> x -> counts_to_histogram(x, bin_edges).weights

    # outside of the fit and in the blinded region the bins are empty
    @test hists[:total].edges[1][1] == 100.0
    @test hists[:total].edges[1][end] == 900.0
    @test any(==(0), hists[:total].weights)
    @test sum(hists[:total].weights) ≈ sum(model(parameters))

    # the order of the parameters does not matter
    reversed = get_best_fit_histograms((b = 4.0, a = 2.0), model, bin_edges)
    @test reversed[:total].weights ≈ hists[:total].weights

    @test_throws ArgumentError get_best_fit_histograms((total = 1.0,), model, bin_edges)

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
    # a flat distribution from 0 to 100 and one from 100 to 200, with weights 1 and 3
    samples = (
        v = [(a = Float64(x), b = 100.0 + x) for x = 0:100],
        weight = [x < 50 ? 1 : 3 for x = 0:100],
    )

    intervals = get_intervals(samples)
    @test Set(keys(intervals)) == Set([:a, :b])
    low, high = intervals[:a]
    @test low < high
    @test intervals[:b] == (100 + intervals[:a][1], 100 + intervals[:a][2])

    # the weights move the interval to the heavier values
    unweighted = get_intervals((v = samples.v, weight = ones(Int, length(samples.v))))
    @test unweighted[:a][1] < intervals[:a][1]

    # a wider level gives a wider interval
    wide = get_intervals(samples; level = 0.95)
    @test wide[:a][1] < low && wide[:a][2] > high
    @test_throws ArgumentError get_intervals(samples; level = 1.0)

    file = tempname() * ".yaml"
    parameters = (a = 60.0, b = 160.0)
    write_activities(file, parameters, samples)
    result = YAML.load_file(file)
    @test Set(keys(result)) == Set(["a", "b"])
    @test result["a"]["best_fit"] == 60.0
    @test result["a"]["low"] ≈ low
    @test result["a"]["high"] ≈ high
end
