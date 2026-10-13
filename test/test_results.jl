using Refit
using StatsBase
using Test

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
    @test hists[:a].weights ≈
          model((a = 2.0, b = 0.0)) |> x -> counts_to_histogram(x, bin_edges).weights

    # outside of the fit and in the blinded region the bins are empty
    @test hists[:total].edges[1][1] == 100.0
    @test hists[:total].edges[1][end] == 900.0
    @test any(==(0), hists[:total].weights)
    @test sum(hists[:total].weights) ≈ sum(model(parameters))

    # the order of the parameters does not matter
    reversed = get_best_fit_histograms((b = 4.0, a = 2.0), model, bin_edges)
    @test reversed[:total].weights ≈ hists[:total].weights

    @test_throws ArgumentError get_best_fit_histograms((total = 1.0,), model, bin_edges)
end
