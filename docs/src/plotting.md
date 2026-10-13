# Plotting

Plots of the results are made with two packages. The spectra, the best fit and its
sources are drawn with [CairoMakie](https://docs.makie.org/), by the functions of this
package. The posterior distributions are drawn with the plot recipes of
[BAT.jl](https://bat.github.io/), which use [Plots.jl](https://docs.juliaplots.org/).
The [tutorial](tutorial.md) shows both for a complete fit.

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

## Histograms of the best fit

The plots compare the measured spectrum with the prediction of the best fit as
histograms in the bins of the fit. `get_best_fit_histograms` makes the histograms of
the prediction for the parameters `best` (for example `mode(samples)`). The result has
the total in the key `:total`, followed by the contribution of each source. The bins
outside of the fit range and in the blinded region are empty.
`counts_to_histogram` turns the measured counts, in the same bins, into a histogram:

```@example demo
data_hist = counts_to_histogram(get_counts(data, bin_edges), bin_edges)
hists = get_best_fit_histograms(best, model, bin_edges)

keys(hists)
```

## The plot of the best fit

`plot_best_fit` draws the measured spectrum as a filled histogram, the total of the best
fit as a black line and each source as a line with its own color. Below it are the
residuals in units of standard deviations ([`normalised_poisson_residual`](@ref)), and
the empty bins are a grey band. The result is a CairoMakie `Figure`.

```@example demo
fig = plot_best_fit(data_hist, hists; title = "Best fit")
```

It can be saved with `CairoMakie.save("fit.png", fig)`. `scale = :linear` gives a
linear y axis instead of the logarithmic one.

To plot fewer lines, `group_histograms` sums sources into groups. A group is a label
with the sources that go into it: a list of names, a regular expression that is matched
to the names, or `:rest` for all the others. The result has the same format as `hists`:

```@example demo
grouped = group_histograms(
    hists,
    ["K40" => r"K40", "Pb212" => [:Pb212_to_Pb208_hpge_bulk], "other" => :rest],
)

plot_best_fit(data_hist, grouped)
```

## Plots with a different style

`plot_best_fit` is made of a few functions that draw on a CairoMakie axis, and they can
be used to make other plots:

- `fit_axes` makes the figure, with an axis for the spectra and an axis for the
  residuals below it. `scale`, `title` and `ylabel` set the y axis and the title.
  It returns `fig, ax, ax_res`.
- `plot_hist!(ax, h)` draws the histogram `h` as a line, or as a filled area with
  `filled = true`. The other keywords, for example `color`, `linewidth`, `linestyle`
  and `label`, go to the CairoMakie function that draws it. Without `color` the lines
  get the next color of the palette.
- `mark_blinded!(ax, h)` draws the grey band over the empty bins of `h`.
- `plot_residuals!(ax, prediction, data)` draws the residuals as points, with the bands
  of 1, 2 and 3 standard deviations unless `bands = false`.

For example, the data as a black line, the best fit in red and the sum of the other
sources dashed:

```@example demo
fig, ax, ax_res = fit_axes(; ylabel = "Counts per 10 keV")

plot_hist!(ax, data_hist; color = :black, label = "Data")
plot_hist!(ax, hists[:total]; color = :red, linewidth = 3, label = "Best fit")
plot_hist!(ax, hists[:K40_hpge_bulk]; linestyle = :dash, label = "K40")

mark_blinded!(ax, hists[:total])
mark_blinded!(ax_res, hists[:total])
plot_residuals!(ax_res, hists[:total], data_hist; bands = false, color = :red)

CairoMakie.xlims!(ax, 500, 3200)
CairoMakie.ylims!(ax_res, -6, 6)
CairoMakie.axislegend(ax)
fig
```

These functions call CairoMakie only with its full names, for example
`CairoMakie.lines!`, and this package does not export any CairoMakie function. If you load
both CairoMakie and Plots.jl with `using`, names such as `plot` and `scatter` clash, so
load one of them with `import`.

## Posterior distributions

The samples of the posterior from `bat_sample` are plotted with the plot recipes of
BAT.jl, which are used with Plots.jl. Load `Plots` and call `plot` on the samples. The
names of the parameters are the names of the activities:

```julia
using BAT, Plots

# the marginal distribution of one parameter, with the mode as a dashed line
plot(samples, :cables_K40; globalmode = true, nbins = 50)

# the distribution of a pair of parameters, to see their correlation
plot(samples, (:cables_K40, :fibers_K40); nbins = 50)

# the distributions of all parameters and of all pairs
plot(samples; nbins = 50)
```

The colored regions are the smallest intervals that contain about 68%, 95% and 99.7% of
the distribution. The keywords `mean`, `std`, `globalmode` and `marginalmode` mark the
mean with its standard deviation, the global mode and the mode of the marginal
distribution. `xlabel`, `ylabel` and `title` are the usual keywords of Plots.jl. See the
[plotting page of BAT.jl](https://bat.github.io/BAT.jl/stable/plotting/) for all the
options. BAT.jl also computes `mode(samples)`, `mean(samples)` and `std(samples)`, and
the smallest credible intervals, see [Saving the results](saving.md).
