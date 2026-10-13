import CairoMakie

# step coordinates of a histogram: the weight of each bin is repeated at both edges
function _steps(h::Histogram{<:Real,1}, floor_value)
    edges = collect(h.edges[1])
    x = repeat(edges; inner = 2)[2:(end-1)]
    y = repeat(max.(Float64.(h.weights), floor_value); inner = 2)
    return x, y
end

# the lowest value that is drawn: a logarithmic axis cannot show 0
_floor_value(ax) = ax.yscale[] === log10 ? 1e-9 : 0.0

"""
    plot_hist!(ax, h; filled = false, kwargs...)

Draw the `Histogram` `h` as a step histogram on the CairoMakie axis `ax`: a line, or a
filled area with `filled = true`. The other keyword arguments, for example `color`,
`linewidth` or `label`, go to the CairoMakie function that draws it (`lines!` or
`band!`). Empty bins are drawn at 0, or at a very small value on a logarithmic axis.
Returns the plot object.
"""
function plot_hist!(ax, h::Histogram{<:Real,1}; filled::Bool = false, kwargs...)
    floor_value = _floor_value(ax)
    x, y = _steps(h, floor_value)
    filled && return CairoMakie.band!(ax, x, fill(floor_value, length(y)), y; kwargs...)
    return CairoMakie.lines!(ax, x, y; kwargs...)
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
    attributes = merge((color = :black, markersize = 6), NamedTuple(kwargs))
    return CairoMakie.scatter!(
        ax,
        centers[fitted],
        normalised_poisson_residual(prediction.weights[fitted], data.weights[fitted]);
        attributes...,
    )
end

"""
    fit_axes(; scale = :log, title = "", ylabel = "Counts")

Make a CairoMakie figure for the plot of a fit: an axis for the spectra with an axis for
the residuals below it. The energy axes are linked and only the lower one has labels.
Returns `(fig, ax, ax_residuals)`. With `scale = :log` the y axis of `ax` is
logarithmic, with `:linear` it is linear. Draw on the axes with [`plot_hist!`](@ref),
[`mark_blinded!`](@ref) and [`plot_residuals!`](@ref), and set the limits, for example
with `CairoMakie.xlims!(ax, 100, 4000)`.
"""
function fit_axes(;
    scale::Symbol = :log,
    title::AbstractString = "",
    ylabel::AbstractString = "Counts",
)
    scale in (:log, :linear) || throw(ArgumentError("scale must be :log or :linear"))

    fig = CairoMakie.Figure(size = (1300, 600))
    ax =
        CairoMakie.Axis(fig[1, 1]; ylabel, yscale = scale == :log ? log10 : identity, title)
    ax_res = CairoMakie.Axis(fig[2, 1]; ylabel = "Residual [σ]", xlabel = "Energy [keV]")

    CairoMakie.rowsize!(fig.layout, 1, CairoMakie.Relative(0.7))
    CairoMakie.rowgap!(fig.layout, 5)
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
the groups, which are named by their label as a `Symbol`. It can be used in place of
`hists`, for example in [`plot_best_fit`](@ref). A source goes into the first group that
selects it, and a group without sources is left out. A name that does not exist, or a
regular expression that matches no source, is an error.

```julia
grouped = group_histograms(hists, ["K40" => r"K40", "Pb212" => [:Pb212_cables], "other" => :rest])
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

"""
    plot_best_fit(data, hists; scale = :log, title = "")

Plot the measured spectrum, the best fit and its sources, with the residuals below.
Returns a CairoMakie `Figure`, which can be saved with `CairoMakie.save("fit.png", fig)`.

`data` is the measured spectrum as a `Histogram` (see [`counts_to_histogram`](@ref))
and `hists` is the result of [`get_best_fit_histograms`](@ref), or of
[`group_histograms`](@ref) to plot fewer lines, with the same bins. The data are the
filled histogram, `hists[:total]` the black line and each of the other histograms a
line with its own color, named in the legend by its key. The empty bins of
`hists[:total]`, the blinded region, are a grey band and have no residual. The energy
axis covers the range of `hists[:total]`. With `scale = :log` the y axis is logarithmic,
and with `:linear` it is linear.

The plot is made of [`fit_axes`](@ref), [`plot_hist!`](@ref), [`mark_blinded!`](@ref)
and [`plot_residuals!`](@ref), which can be used to make a plot with a different style.
"""
function plot_best_fit(
    data::Histogram{<:Real,1},
    hists::AbstractDict{Symbol,<:Histogram{<:Real,1}};
    scale::Symbol = :log,
    title::AbstractString = "",
)
    haskey(hists, :total) ||
        throw(ArgumentError("hists must have the key :total, see get_best_fit_histograms"))
    total = hists[:total]
    data.edges == total.edges || throw(ArgumentError("data and hists differ in binning"))

    bin_width = round(minimum(diff(total.edges[1])), sigdigits = 3)
    fig, ax, ax_res = fit_axes(; scale, title, ylabel = "Counts per $bin_width keV")

    # the data, the sources (with the colors of Makie's palette) and the best fit
    plot_hist!(ax, data; filled = true, color = :lightblue, label = "Data")
    for (name, h) in hists
        name == :total || plot_hist!(ax, h; label = string(name))
    end
    plot_hist!(ax, total; color = :black, linewidth = 2, label = "Best fit")

    mark_blinded!(ax, total)
    mark_blinded!(ax_res, total)
    plot_residuals!(ax_res, total, data)

    edges = total.edges[1]
    CairoMakie.xlims!(ax, first(edges), last(edges))
    CairoMakie.ylims!(
        ax,
        scale == :log ? 0.5 : 0,
        (scale == :log ? 3 : 1.1) * maximum(data.weights),
    )
    CairoMakie.ylims!(ax_res, -6, 6)
    CairoMakie.axislegend(ax; nbanks = 3)
    CairoMakie.axislegend(ax_res; nbanks = 3)

    return fig
end
