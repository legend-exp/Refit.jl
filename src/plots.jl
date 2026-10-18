import CairoMakie

# step coordinates of a histogram: the weight of each bin is repeated at both edges,
# and the line starts and ends at `floor_value`, so that the first and the last bin
# are closed
function _steps(h::Histogram{<:Real,1}, floor_value)
    edges = collect(h.edges[1])
    x = repeat(edges; inner = 2)
    y = [
        floor_value;
        repeat(max.(Float64.(h.weights), floor_value); inner = 2);
        floor_value
    ]
    return x, y
end

# the lowest value that is drawn: a logarithmic axis cannot show 0
_floor_value(ax) = ax.yscale[] === log10 ? 1e-9 : 0.0

"""
    plot_hist!(ax, h; filled = false, kwargs...)

Draw the `Histogram` `h` as a step histogram on the CairoMakie axis `ax`: a line (with
`linewidth = 2` unless set), or a filled area with `filled = true`. The other keyword
arguments, for example `color`, `linewidth` or `label`, go to the CairoMakie function
that draws it (`lines!` or `band!`). Empty bins are drawn at 0, or at a very small value
on a logarithmic axis. Returns the plot object.

With bins of equal width the counts per bin are drawn. With bins of different widths,
for example from [`uniform_binning_with_peaks`](@ref) or a blinded region, each bin is
divided by its width, so the plot shows counts per unit of the x axis (counts per keV
for a spectrum in keV).
"""
function plot_hist!(ax, h::Histogram{<:Real,1}; filled::Bool = false, kwargs...)
    widths = diff(h.edges[1])
    if !all(w -> w ≈ first(widths), widths)
        h = Histogram(h.edges, h.weights ./ widths)
    end
    floor_value = _floor_value(ax)
    x, y = _steps(h, floor_value)
    filled && return CairoMakie.band!(ax, x, fill(floor_value, length(y)), y; kwargs...)
    attributes = merge((linewidth = 2,), NamedTuple(kwargs))
    return CairoMakie.lines!(ax, x, y; attributes...)
end

"""
    mark_blinded!(ax, h; kwargs...)

Mark the empty bins of the `Histogram` `h` on the CairoMakie axis `ax` with a vertical
band. These are the bins that are not fitted, for example the blinded region in the
histograms of [`get_best_fit_histograms`](@ref). The keyword arguments, for example
`color`, go to `CairoMakie.vspan!`, and the default is a grey band.
"""
function mark_blinded!(ax, h::Histogram{<:Real,1}; kwargs...)
    attributes = merge((color = (:grey, 0.5),), NamedTuple(kwargs))
    edges = h.edges[1]
    for i in findall(==(0), h.weights)
        CairoMakie.vspan!(ax, edges[i], edges[i+1]; attributes...)
    end
    return nothing
end

"""
    plot_residuals!(ax, prediction, data; bands = true, kwargs...)

Draw the residuals of the `data` with respect to the `prediction` (both `Histogram`s
with the same bins) as points on the CairoMakie axis `ax`, in units of standard
deviations, see [`normalised_poisson_residual`](@ref). Bins where the prediction is
empty are not fitted and are left out. With `bands = true` the bands of 1, 2 and 3
standard deviations are drawn behind the points. The keyword arguments, for example
`color` or `markersize`, go to `CairoMakie.scatter!`. Returns the plot object.
"""
function plot_residuals!(
    ax,
    prediction::Histogram{<:Real,1},
    data::Histogram{<:Real,1};
    bands::Bool = true,
    kwargs...,
)
    prediction.edges == data.edges ||
        throw(ArgumentError("prediction and data differ in binning"))

    if bands
        CairoMakie.hspan!(ax, -3, 3; color = (:lightcyan, 0.6), label = "3σ")
        CairoMakie.hspan!(ax, -2, 2; color = (:lightskyblue, 0.3), label = "2σ")
        CairoMakie.hspan!(ax, -1, 1; color = (:dodgerblue2, 0.5), label = "1σ")
    end

    fitted = findall(!=(0), prediction.weights)
    edges = prediction.edges[1]
    centers = (edges[1:(end-1)] .+ edges[2:end]) ./ 2
    attributes =
        merge((color = :black, marker = :circle, markersize = 10), NamedTuple(kwargs))
    return CairoMakie.scatter!(
        ax,
        centers[fitted],
        normalised_poisson_residual(prediction.weights[fitted], data.weights[fitted]);
        attributes...,
    )
end

# The shortest interval that contains the fraction `level` of the weights of the samples `x`.
function _smallest_interval(
    x::AbstractVector{<:Real},
    weights::AbstractVector{<:Real},
    level::Real,
)
    length(x) == length(weights) ||
        throw(DimensionMismatch("samples and weights differ in length"))
    order = sortperm(x)
    sorted = x[order]
    cumulative = cumsum(Float64.(weights[order]))
    target = level * last(cumulative)

    best_width = Inf
    best = (first(sorted), last(sorted))
    stop = 1
    for start in eachindex(sorted)
        # the weight in start:stop is cumulative[stop] - cumulative[start - 1]
        before = start == 1 ? 0.0 : cumulative[start-1]
        while cumulative[stop] - before < target
            stop == length(sorted) && return best
            stop += 1
        end
        width = sorted[stop] - sorted[start]
        if width < best_width
            best_width = width
            best = (sorted[start], sorted[stop])
        end
    end
    return best
end

"""
    plot_posterior!(ax, samples, name; level = 0.68, nbins = 50, best = nothing,
                    truth = nothing, kwargs...)

Draw the marginal posterior distribution of the parameter `name` on the CairoMakie axis
`ax`, as a step histogram normalised to an integral of 1. `samples` are the samples of the
posterior, for example from `bat_sample`: the parameters of each sample are in
`samples.v` and its weight in `samples.weight`. The histogram has `nbins` bins between the
smallest and the largest value.

The bins inside the smallest interval that contains the fraction `level` of the
distribution are filled. This is the same interval as the one saved by
[`write_parameter_values`](@ref) when it comes from `BAT.smallest_credible_intervals`.
The filled region has the color of the histogram with some transparency, or `bandcolor`,
for example `bandcolor = :lightblue`. With `best` (for example `mode(samples)[name]`) a
dashed vertical line marks the best fit, and with `truth` a red line marks the true
value. The lower limit of the y axis is 0. The other keywords, for example `color` or
`label`, go to [`plot_hist!`](@ref).

The plots have labels, so a legend can be added with `CairoMakie.axislegend(ax)`, or in
another place with `CairoMakie.Legend(fig[2, 4], ax)`.

Returns the interval as `low .. high`.
"""
function plot_posterior!(
    ax,
    samples,
    name::Symbol;
    level::Real = 0.68,
    nbins::Integer = 50,
    best::Union{Nothing,Real} = nothing,
    truth::Union{Nothing,Real} = nothing,
    bandcolor = nothing,
    kwargs...,
)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    haskey(first(samples.v), name) || throw(ArgumentError("unknown parameter $name"))

    values = [v[name] for v in samples.v]
    weights = Float64.(samples.weight)
    # the last bin includes the largest value
    lowest, highest = extrema(values)
    lowest < highest || throw(ArgumentError("all the samples of $name have the same value"))
    edges = range(lowest, highest + 1e-9 * (highest - lowest); length = nbins + 1)

    histogram = fit(Histogram, values, Weights(weights), edges)
    density = histogram.weights ./ (sum(histogram.weights) * step(edges))
    posterior = Histogram(edges, density)

    low, high = _smallest_interval(values, weights, level)
    centers = (edges[1:(end-1)] .+ edges[2:end]) ./ 2
    inside = (low .<= centers) .& (centers .<= high)

    attributes = merge((color = :black, label = "Posterior"), NamedTuple(kwargs))
    fill_color = isnothing(bandcolor) ? (attributes.color, 0.4) : bandcolor
    plot_hist!(
        ax,
        Histogram(edges, density .* inside);
        filled = true,
        color = fill_color,
        label = "Smallest $(round(Int, 100 * level))% interval",
    )
    plot_hist!(ax, posterior; attributes...)
    isnothing(best) || CairoMakie.vlines!(
        ax,
        best;
        color = attributes.color,
        linestyle = :dash,
        label = "Best fit",
    )
    isnothing(truth) || CairoMakie.vlines!(ax, truth; color = :red, label = "True value")
    CairoMakie.ylims!(ax; low = 0)  # a density cannot be negative

    return low .. high
end

"""
    fit_axes(; scale = :log, title = "", ylabel = "Counts", size = (1300, 600), colors)

Make a CairoMakie figure for the plot of a fit: an axis for the spectra with an axis for
the residuals below it. The energy axes are linked and only the lower one has labels.
Returns `(fig, ax, ax_residuals)`.

`size` is the size of the figure in pixels. With `scale = :log` the y axis of `ax` is
logarithmic, with `:linear` it is linear. The lines drawn by [`plot_hist!`](@ref)
without a `color` take the next of `colors`, a vector of colors that repeats after the
last one. The default is a color-blind friendly palette of 7 colors.

Draw on the axes with [`plot_hist!`](@ref), [`mark_blinded!`](@ref) and
[`plot_residuals!`](@ref), and set the limits, for example with
`CairoMakie.xlims!(ax, 100, 4000)`.
"""
function fit_axes(;
    scale::Symbol = :log,
    title::AbstractString = "",
    ylabel::AbstractString = "Counts",
    size::Tuple{Integer,Integer} = (1300, 600),
    colors::AbstractVector = CairoMakie.Makie.wong_colors(),
)
    scale in (:log, :linear) || throw(ArgumentError("scale must be :log or :linear"))

    fig = CairoMakie.Figure(; size)
    ax = CairoMakie.Axis(
        fig[1, 1];
        ylabel,
        yscale = scale == :log ? log10 : identity,
        title,
        palette = (color = colors,),
    )
    ax_res = CairoMakie.Axis(fig[2, 1]; ylabel = "Residual [σ]", xlabel = "Energy [keV]")

    CairoMakie.rowsize!(fig.layout, 1, CairoMakie.Relative(0.7))
    CairoMakie.rowgap!(fig.layout, 12)
    CairoMakie.hidexdecorations!(ax; grid = false)
    CairoMakie.linkxaxes!(ax, ax_res)
    return fig, ax, ax_res
end

"""
    group_histograms(hists, groups)

Sum the sources of the best fit into groups, for example to plot fewer lines. `hists` is
the result of [`get_best_fit_histograms`](@ref) and `groups` is a vector with one
`label => selection` for each group. The selection is one of
- a vector of source names, for example `[:K40_cables, :K40_fibers]`,
- a regular expression that is matched to the source names, for example `r"K40"`,
- `:rest`, for all the sources that no other group selects.

Returns an ordered dictionary in the same format as `hists`, with `:total` followed by
the groups, which are named by their label as a `Symbol`. The label can be any text,
including Unicode such as `"²¹²Pb"` or `"2νββ"`. A source goes into the first
group that selects it, and a group without sources is left out. A name that does not
exist, or a regular expression that matches no source, is an error.

```julia
grouped = group_histograms(
    hists,
    ["K40" => r"K40", "Pb212" => [:Pb212_cables], "other" => :rest],
)
```
"""
function group_histograms(
    hists::AbstractDict{Symbol,<:Histogram{<:Real,1}},
    groups::AbstractVector{<:Pair},
)
    haskey(hists, :total) ||
        throw(ArgumentError("hists must have the key :total, see get_best_fit_histograms"))
    sources = [name for name in keys(hists) if name != :total]

    selected = Dict{Symbol,Vector{Symbol}}()
    taken = Set{Symbol}()
    for (label, selection) in groups
        selection == :rest && continue
        if selection isa Regex
            members = filter(s -> occursin(selection, string(s)), sources)
            isempty(members) &&
                throw(ArgumentError("no source matches $selection for group \"$label\""))
        elseif selection isa AbstractVector && all(s -> s isa Symbol, selection)
            members = collect(Symbol, selection)
            unknown = setdiff(members, sources)
            isempty(unknown) || throw(ArgumentError("unknown sources in groups: $unknown"))
        else
            throw(
                ArgumentError(
                    "the selection of group \"$label\" must be a vector of names, a Regex or :rest",
                ),
            )
        end
        label isa Union{Symbol,AbstractString} ||
            throw(ArgumentError("the label of a group must be a String or a Symbol"))
        selected[Symbol(label)] = filter(!in(taken), members)
        union!(taken, members)
    end

    grouped = OrderedDict{Symbol,Histogram{Float64,1}}(:total => hists[:total])
    for (label, selection) in groups
        name = Symbol(label)
        name == :total && throw(ArgumentError("a group cannot be called total"))
        members = selection == :rest ? setdiff(sources, taken) : selected[name]
        isempty(members) || (grouped[name] = add_histograms([hists[s] for s in members]))
    end
    return grouped
end
