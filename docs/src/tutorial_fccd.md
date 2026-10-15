# Tutorial: fitting the FCCD with template morphing

This tutorial fits the dead layer of a detector to the spectrum of a 60 keV γ line, with
spectra simulated on a grid of dead layer models (see
[Accounting for systematic uncertainties](generalised_hist.md)). The data are made with known values, so
the result can be compared with the truth at the end.

## Model of the detector surface

A γ of 60 keV deposits its energy at a depth that is uniform from 0 to 2 mm below the
surface of the detector. The fraction of the charge that is collected from a depth ``d``
depends on two parameters, in mm:

- the FCCD (full charge collection depth): below it, all the charge is collected,
- the thickness `TL` of the transition layer above the FCCD: through it, the collected
  fraction rises linearly from 0 at a depth of FCCD − TL to 1 at the FCCD. Above it is
  the dead layer, from which no charge is collected.

The measured energy is 60 keV times the collected fraction. So the spectrum has two
parts:

- a peak at 60 keV, from depths below the FCCD, with probability (2 mm − FCCD) / 2 mm,
- a flat continuum from 0 to 60 keV, from the transition layer, with probability
  TL / 2 mm.

Events from the dead layer have no energy and are not recorded. The counts in each bin
are linear in the FCCD and in the TL, so here the interpolation between grid points is
exact, apart from the statistical fluctuations of the simulation.

## Simulated templates

The input files are in `test/test_files/fccd` of the repository. The Python script
`make_files.py` in the same directory makes them, with the same tools that legend-simflow
uses for its PDFs.

- `templates.lh5` has the spectra of 10⁶ simulated γs at each point of a grid of 5 FCCD
  values from 0.6 to 1.4 mm, times 6 TL values from 0 to 0.5 mm, in 1 keV bins. They form
  one generalised histogram, `pdf/geds`. The file also holds the number of simulated γs.
- `data.lh5` has the measured spectrum, `hist/geds`.

```@example fccd
dir = "../../test/test_files/fccd"  # relative to docs/src
readdir(dir)
```

`read_generalised_pdf` reads the templates and divides them by the number of simulated
γs, so they are in counts per γ. The result knows the grid, in mm:

```@example fccd
using Refit
using StatsBase

templates = read_generalised_pdf(joinpath(dir, "templates.lh5"); groups = "pdf/geds")
templates.grid
```

`get_histogram` gives the spectrum at any point inside the grid. The peak falls and the
continuum grows as the dead layer gets thicker:

```@example fccd
using LegendMakie
import Format
import CairoMakie

fig = CairoMakie.Figure(size = (900, 450))
ax = CairoMakie.Axis(
    fig[1, 1];
    yscale = log10,
    xlabel = "Energy [keV]",
    ylabel = "Counts per γ and keV",
)
for (fccd, tl) in ((0.6, 0.1), (1.0, 0.3), (1.4, 0.5))
    h = get_histogram(templates, (fccd = fccd, tl = tl))
    plot_hist!(ax, h; label = "FCCD = $fccd mm, TL = $tl mm")
end
CairoMakie.ylims!(ax, 1e-4, 1)
CairoMakie.axislegend(ax; position = :lt)
fig
```

## Binning

The fit uses 5 keV bins from 5 to 65 keV, which leaves out the lowest bins, as an
energy threshold does. Rebinning to these edges keeps only this range, for the
templates and for the data alike:

```@example fccd
edges = collect(5.0:5.0:65.0)
templates_fit = rebin_hist(templates, edges)

data = rebin_hist(read_hist("hist/geds", joinpath(dir, "data.lh5")), edges)
counts = data.weights
nothing # hide
```

## Model and likelihood

The source emitted 20000 γs during the measurement. This number must be known, for
example from the activity of the source and the livetime. If it were a free
parameter, the fit could not tell a thicker dead layer from fewer γs: the spectrum only
measures the height of the peak and of the continuum.

The model gives the expected counts in each bin: the number of γs times the counts per
γ, interpolated at the FCCD and TL of the parameters. `build_likelihood` takes any such
function of the parameters.

```@example fccd
n_gamma = 20_000

model(p) = n_gamma .* get_weights(templates_fit, p)
likelihood = build_likelihood(counts, model)
nothing # hide
```

## Priors and sampling

The prior of each parameter is uniform over its grid, the range where the templates can
be interpolated.

```@example fccd
using IntervalSets
using BAT

prior = build_prior_dist(Dict(:fccd => 0.6 .. 1.4, :tl => 0.0 .. 0.5))

samples = bat_sample(
    PosteriorMeasure(likelihood, prior),
    TransformedMCMC(
        proposal = RandomWalk(),
        strict = false,
        burnin = MCMCMultiCycleBurnin(nsteps_per_cycle = 20_000, max_ncycles = 50),
        nsteps = 10^5,
        nchains = 4,
    ),
).result
nothing # hide
```

## Result

The fake data were made with an FCCD of 1.05 mm and a TL of 0.27 mm. The marginal
posterior of each parameter, with the best fit and the true value:

```@example fccd
truth = (fccd = 1.05, tl = 0.27)  # mm
best = mode(samples)

fig = CairoMakie.Figure(size = (900, 400))
for (i, (name, title)) in enumerate(((:fccd, "FCCD"), (:tl, "Transition layer")))
    ax = CairoMakie.Axis(
        fig[1, i];
        title,
        xlabel = "Thickness [mm]",
        ylabel = "Probability density",
        xticks = CairoMakie.LinearTicks(3),
    )
    plot_posterior!(ax, samples, name; best = best[name], truth = truth[name])
end
CairoMakie.Legend(fig[1, 3], fig.content[1]; framevisible = false)
fig
```

The prediction of the best fit compared with the data:

```@example fccd
bin_edges = (low = edges[1:(end-1)], high = edges[2:end])
best_hist = counts_to_histogram(model(best), bin_edges)

fig, ax, ax_res = fit_axes(; ylabel = "Counts per 5 keV")
plot_hist!(ax, data; filled = true, alpha = 0.6, color = :lightblue, label = "Data")
plot_hist!(ax, best_hist; color = :black, linewidth = 3, label = "Best fit")
plot_residuals!(ax_res, best_hist, data)
CairoMakie.ylims!(ax, 100, 1e5)
CairoMakie.axislegend(ax; position = :lt)
fig
```
