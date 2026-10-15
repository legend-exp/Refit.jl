import CairoMakie
using Distributions
using IntervalSets
using Random
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
        # 400 to 600 keV left out, as a blinded region
        binning = [100:100:400, 600:100:900],
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

    clash, clash_edges = build_forward_model(Dict(:total => pdfs[:a]); e_low = 100)
    @test_throws ArgumentError get_best_fit_histograms((total = 1.0,), clash, clash_edges)

    # components with an activity function and a shared parameter
    shared, shared_edges = build_forward_model(
        Dict(
            :a => (pdf = pdfs[:a], activity = p -> p.s),
            :b => (pdf = pdfs[:b], activity = p -> 2 * p.s),
        );
        e_low = 100,
        e_high = 900,
    )
    parts = get_best_fit_histograms((s = 2.0,), shared, shared_edges)
    @test Set(keys(parts)) == Set([:total, :a, :b])
    @test parts[:total].weights ≈ parts[:a].weights .+ parts[:b].weights

    # an analytic component is drawn like the others
    with_flat, flat_edges = build_forward_model(
        Dict(:a => pdfs[:a], :flat => (analytic = (E, p) -> p.rate,));
        e_low = 100,
        e_high = 900,
    )
    flat_parts = get_best_fit_histograms((a = 1.0, rate = 0.5), with_flat, flat_edges)
    @test Set(keys(flat_parts)) == Set([:total, :a, :flat])
    @test flat_parts[:total].weights ≈ flat_parts[:a].weights .+ flat_parts[:flat].weights
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

@testset "fit_axes" begin
    fig, ax, ax_res = fit_axes(; title = "test", ylabel = "Counts per 20 keV")
    @test fig isa CairoMakie.Figure
    @test ax isa CairoMakie.Axis
    @test ax_res isa CairoMakie.Axis
    @test ax.yscale[] === log10
    @test ax.title[] == "test"

    # the size of the figure
    fig_big, _, _ = fit_axes(; size = (800, 700))
    @test Tuple(fig_big.scene.viewport[].widths) == (800, 700)

    # lines without a color cycle through the colors
    _, ax_colors, _ = fit_axes(; colors = [:red, :blue])
    h = Histogram(0.0:1.0:3.0, [1.0, 2.0, 3.0])
    @test plot_hist!(ax_colors, h).color[] == CairoMakie.to_color(:red)
    @test plot_hist!(ax_colors, h).color[] == CairoMakie.to_color(:blue)

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
        # 400 to 600 keV left out, as a blinded region
        binning = [100:100:400, 600:100:900],
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

    # the lines have a default width that can be changed
    @test plot_hist!(ax, total).linewidth[] == 2
    @test plot_hist!(ax, total; linewidth = 5).linewidth[] == 5
    @test plot_hist!(ax_log, data; filled = true) isa CairoMakie.Makie.AbstractPlot

    # the step line goes up at the first edge and down at the last one
    x, y = Refit._steps(Histogram((0.0:1.0:2.0,), [3.0, 4.0]), 0.0)
    @test x == [0, 0, 1, 1, 2, 2]
    @test y == [0, 3, 3, 4, 4, 0]

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

@testset "plot_posterior!" begin
    # the shortest interval of a skewed distribution starts at its lowest values
    x = collect(0.0:99.0)
    weights = 100.0 .- x
    low, high = Refit._smallest_interval(x, weights, 0.68)
    @test low == 0.0
    # brute force: the shortest of all the intervals that contain 68% of the weight
    widths = [
        (x[j] - x[i], i, j) for i in eachindex(x) for
        j = i:length(x) if sum(weights[i:j]) >= 0.68 * sum(weights)
    ]
    @test (high - low) == minimum(first.(widths))

    # a symmetric distribution: the same as the central interval, within the bin
    Random.seed!(1)
    gauss = randn(100_000)
    low, high = Refit._smallest_interval(gauss, ones(length(gauss)), 0.6827)
    @test low ≈ -1 atol = 0.05
    @test high ≈ 1 atol = 0.05
    @test_throws DimensionMismatch Refit._smallest_interval([1.0, 2.0], [1.0], 0.68)

    # samples in the format of bat_sample: the parameters in v and the weights
    samples = (v = [(a = x, b = 2x) for x in gauss[1:20_000]], weight = rand(1:3, 20_000))
    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1, 1])
    interval = plot_posterior!(ax, samples, :a; best = 0.0, truth = 0.5, color = :blue)
    @test ax.limits[][2][1] == 0  # the y axis starts at 0
    @test interval isa ClosedInterval
    @test interval.left < 0 < interval.right
    @test interval.right - interval.left ≈ 2 atol = 0.15

    # the other parameter is twice as wide, and a different level and number of bins
    ax_b = CairoMakie.Axis(fig[2, 1])
    interval_b =
        plot_posterior!(ax_b, samples, :b; level = 0.95, nbins = 20, bandcolor = :lightblue)
    @test interval_b.right - interval_b.left ≈ 8 atol = 0.5

    # the filled band has the color of the line with transparency, or bandcolor
    bands(axis) = [p for p in axis.scene.plots if p isa CairoMakie.Makie.Band]
    @test only(bands(ax)).color[] == (:blue, 0.4)
    @test only(bands(ax_b)).color[] == :lightblue

    # the plots have labels for a legend
    @test only(bands(ax)).label[] == "Smallest 68% interval"
    @test only(bands(ax_b)).label[] == "Smallest 95% interval"
    legend = CairoMakie.axislegend(ax)
    @test length(legend.entrygroups[][1][2]) == 4  # posterior, interval, best fit, truth

    file = tempname() * ".png"
    CairoMakie.save(file, fig)
    @test filesize(file) > 0

    @test_throws ArgumentError plot_posterior!(ax, samples, :c)
    @test_throws ArgumentError plot_posterior!(ax, samples, :a; level = 1.0)
    constant = (v = [(a = 1.0,), (a = 1.0,)], weight = [1, 1])
    @test_throws ArgumentError plot_posterior!(ax, constant, :a)
end
