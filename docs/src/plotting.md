# Plotting

The spectra, the best fit and its sources are plotted with
[CairoMakie](https://docs.makie.org/), by the functions of this package. The posterior
distributions are plotted with the plot recipes of [BAT.jl](https://bat.github.io/),
which use [Plots.jl](https://docs.juliaplots.org/). The [tutorial](tutorial.md) shows
both for a complete fit.

```@setup demo
include("example_data.jl")

using Refit
import CairoMakie

# the forward model page, with the true activities as the best fit
norm = 1e-3 * livetime
model, bin_edges = build_forward_model(
    mc_dict;
    e_low = 500,
    e_high = 3200,
    blind = (2015.0, 2065.0),
    norm,
)
best = true_activities
```

Loading [LegendMakie.jl](https://legend-exp.github.io/LegendMakie.jl/stable/) changes
all CairoMakie plots to the LEGEND style. It is not needed, and it also needs `Format`
to be loaded:

```@example demo
using LegendMakie
import Format
```

## The plot of the best fit

`plot_best_fit` draws the measured spectrum, the best fit and its sources, with the
residuals below. It needs two histograms of the same bins. The measured spectrum comes
from `counts_to_histogram`, and `get_best_fit_histograms` makes the histograms of the
prediction for the parameters `best` (for example `mode(samples)`). These are the total
in the key `:total` and the contribution of each source. The bins outside of the fit
range and in the blinded region are empty and shown as a grey band.


```@example demo
data_hist = counts_to_histogram(get_counts(data, bin_edges), bin_edges)
hists = get_best_fit_histograms(best, model, bin_edges)

fig = plot_best_fit(data_hist, hists; title = "Best fit")
```

The result is a CairoMakie `Figure`, which can be saved with
`CairoMakie.save("fit.png", fig)`. To plot fewer lines, `group_histograms` sums sources
into groups. A group is a label with the sources that go into it: a list of names, a
regular expression for the names, or `:rest` for all the others.

```@example demo
grouped = group_histograms(
    hists,
    ["K40" => r"K40", "Pb212" => [:Pb212_to_Pb208_hpge_bulk], "other" => :rest],
)

plot_best_fit(data_hist, grouped)
```

## Plots with a different style

`plot_best_fit` is made of functions that draw on a CairoMakie axis, which can be used to
make other plots. `fit_axes` makes the figure with an axis for the spectra and one for
the residuals. `plot_hist!` draws a histogram as a line, or filled with `filled = true`,
and passes the other keywords (`color`, `linewidth`, `label`, ...) to Makie.
`mark_blinded!` draws the grey band and `plot_residuals!` the residuals. For example:

```@example demo
fig, ax, ax_res = fit_axes(; ylabel = "Counts per 10 keV")

plot_hist!(ax, data_hist; color = :black, label = "Data")
plot_hist!(ax, hists[:total]; color = :red, linewidth = 3, label = "Best fit")
plot_hist!(ax, hists[:K40_hpge_bulk]; linestyle = :dash, label = "K40")

mark_blinded!(ax, hists[:total])
mark_blinded!(ax_res, hists[:total])
plot_residuals!(ax_res, hists[:total], data_hist; bands = false, color = :red)

CairoMakie.axislegend(ax)
fig
```

These functions call CairoMakie only with its full names, so loading both CairoMakie and
Plots.jl with `using` makes names such as `plot` clash. Load one with `import`.

## Posterior distributions

The posterior samples from `bat_sample` are plotted with the plot recipes of BAT.jl,
which use Plots.jl. See the
[plotting section of the BAT.jl manual](https://bat.github.io/BAT.jl/stable/plotting/).
The [tutorial](tutorial.md) shows them for a fit.
