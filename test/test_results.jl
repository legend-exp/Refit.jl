import CairoMakie
using Distributions
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

@testset "normalised_poisson_residual" begin
    # Gaussian regime (mu >= 50): (x - mu) / sqrt(mu)
    @test normalised_poisson_residual(100.0, 100) ≈ 0.0 atol = 1e-10
    @test normalised_poisson_residual(100.0, 110) ≈ 1.0 atol = 1e-10
    @test normalised_poisson_residual(100.0, 90) ≈ -1.0 atol = 1e-10

    # no counts expected and none seen
    @test normalised_poisson_residual(0.3, 0) == 0.0
    @test normalised_poisson_residual(0.3, 0.5) == 0.0

    # low-mu regime: below the mode the residual is negative, above it positive
    @test normalised_poisson_residual(5.0, 1) < 0.0
    @test normalised_poisson_residual(5.0, 10) > 0.0

    # the p-value converted into standard deviations
    let μ = 5.0, x = 10
        prob = 1 - cdf(Poisson(μ), x - 1)
        @test normalised_poisson_residual(μ, x) ≈ quantile(Normal(), 1 - prob) atol = 1e-10
    end
    let μ = 5.0, x = 1
        prob = cdf(Poisson(μ), x)
        @test normalised_poisson_residual(μ, x) ≈ -quantile(Normal(), 1 - prob) atol = 1e-10
    end

    # arrays
    res = normalised_poisson_residual([100.0, 100.0, 5.0], [110, 90, 0])
    @test res isa AbstractVector
    @test res[1] ≈ 1.0 atol = 1e-10
    @test res[2] ≈ -1.0 atol = 1e-10
    @test res[3] < 0.0
    @test normalised_poisson_residual(100.0, 100) isa Float64
    @test normalised_poisson_residual(5, 3) isa Float64
    @test res isa Vector{Float64}
    @test_throws DimensionMismatch normalised_poisson_residual([1.0, 2.0], [1, 2, 3])
    @test_throws MethodError normalised_poisson_residual(1.0, [1, 2])
end

@testset "group_histograms" begin
    edges = 0.0:100.0:1000.0
    pdfs = Dict(
        :K40_cables => Histogram(edges, fill(10.0, 10)),
        :K40_fibers => Histogram(edges, fill(5.0, 10)),
        :alpha => Histogram(edges, fill(1.0, 10)),
    )
    model, bin_edges = build_forward_model(pdfs; e_low = 100, e_high = 900, norm = 1e-3)
    hists = get_best_fit_histograms(
        (K40_cables = 2e4, K40_fibers = 4e4, alpha = 1e3),
        model,
        bin_edges,
    )

    # the best fit histograms are ordered: the total, then the sources in the order of the parameters
    @test collect(keys(hists)) == [:total, :K40_cables, :K40_fibers, :alpha]

    # a regular expression, a list of names and the rest, in a list of mixed types
    grouped = group_histograms(
        hists,
        ["K40" => r"K40", "fibers" => [:K40_fibers], "other" => :rest],
    )
    @test collect(keys(grouped)) == [:total, :K40, :other]  # the K40 sources are all in the first group
    @test grouped[:total] == hists[:total]
    @test grouped[:K40].weights ≈ hists[:K40_cables].weights .+ hists[:K40_fibers].weights
    @test grouped[:other].weights == hists[:alpha].weights

    # the order of the groups is kept, and a group without sources is left out
    grouped = group_histograms(hists, [:other => :rest, "alpha" => [:alpha]])
    @test collect(keys(grouped)) == [:total, :other, :alpha]
    @test grouped[:other].weights ≈ hists[:K40_cables].weights .+ hists[:K40_fibers].weights
    grouped = group_histograms(hists, ["all" => r".", "other" => :rest])
    @test collect(keys(grouped)) == [:total, :all]

    @test_throws ArgumentError group_histograms(hists, ["x" => [:nothing]])
    @test_throws ArgumentError group_histograms(hists, ["x" => r"nothing"])
    @test_throws ArgumentError group_histograms(hists, ["x" => 1.5])
    @test_throws ArgumentError group_histograms(hists, ["total" => r"K40"])
    @test_throws ArgumentError group_histograms(Dict(:a => hists[:alpha]), ["x" => :rest])
end

@testset "plot_best_fit" begin
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
    hists = get_best_fit_histograms((a = 2e4, b = 4e4), model, bin_edges)
    data = counts_to_histogram(round.(model((a = 2e4, b = 4e4))), bin_edges)

    fig = plot_best_fit(data, hists)
    @test fig isa CairoMakie.Figure
    file = tempname() * ".png"
    CairoMakie.save(file, fig)
    @test filesize(file) > 0

    # grouped histograms, a linear scale and a title
    grouped = group_histograms(hists, ["first" => r"^a", "other" => :rest])
    @test plot_best_fit(data, grouped; scale = :linear, title = "test") isa
          CairoMakie.Figure

    @test_throws ArgumentError plot_best_fit(data, Dict(:a => hists[:a]))
    other_bins = counts_to_histogram(fill(1.0, 3), (low = [0.0, 1, 2], high = [1.0, 2, 3]))
    @test_throws ArgumentError plot_best_fit(other_bins, hists)
end

@testset "fit_axes" begin
    fig, ax, ax_res = fit_axes(; title = "test", ylabel = "Counts per 20 keV")
    @test fig isa CairoMakie.Figure
    @test ax isa CairoMakie.Axis
    @test ax_res isa CairoMakie.Axis
    @test ax.yscale[] === log10
    @test ax.title[] == "test"

    _, ax_linear, _ = fit_axes(; scale = :linear)
    @test ax_linear.yscale[] === identity
    @test_throws ArgumentError fit_axes(; scale = :sqrt)
end

@testset "plot blocks" begin
    edges = 0.0:100.0:1000.0
    pdfs = Dict(:a => Histogram(edges, fill(10.0, 10)))
    model, bin_edges = build_forward_model(
        pdfs;
        e_low = 100,
        e_high = 900,
        blind = (450.0, 550.0),
        norm = 1e-3,
    )
    total = counts_to_histogram(model((a = 2e4,)), bin_edges)
    data = counts_to_histogram(round.(model((a = 2e4,))), bin_edges)

    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    ax_log = CairoMakie.Axis(fig[2, 1]; yscale = log10)
    ax_res = CairoMakie.Axis(fig[3, 1])

    # a line or a filled area, on a linear and on a logarithmic axis (with empty bins)
    @test plot_hist!(ax, total; color = :black, linewidth = 2, label = "line") isa
          CairoMakie.Makie.AbstractPlot
    @test plot_hist!(ax, data; filled = true, color = :lightblue) isa
          CairoMakie.Makie.AbstractPlot
    @test plot_hist!(ax_log, total) isa CairoMakie.Makie.AbstractPlot
    @test plot_hist!(ax_log, data; filled = true) isa CairoMakie.Makie.AbstractPlot

    # the empty bins of the prediction are not drawn as residuals
    @test mark_blinded!(ax, total; color = :red) === nothing
    scatter = plot_residuals!(ax_res, total, data; markersize = 10)
    @test length(scatter[1][]) == count(!=(0), total.weights)
    @test plot_residuals!(ax_res, total, data; bands = false) isa
          CairoMakie.Makie.AbstractPlot

    file = tempname() * ".png"
    CairoMakie.save(file, fig)
    @test filesize(file) > 0

    other_bins = counts_to_histogram(fill(1.0, 3), (low = [0.0, 1, 2], high = [1.0, 2, 3]))
    @test_throws ArgumentError plot_residuals!(ax_res, total, other_bins)
end
