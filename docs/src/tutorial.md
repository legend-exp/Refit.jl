# Tutorial: fitting a background spectrum

This tutorial fits a measured energy spectrum with a sum of simulated spectra, from
reading the file, sampling the posterior and plotting the results.
This runs every time the documentation is built.

The fit uses automatically generated fake data (pseudo-experiments). The measured
spectrum is generated from known activities, so the result of the fit can be compared
with the truth at the end.

## Fake data

The input files are in `test/test_files/tutorial` of the repository, in LH5 format. They
stand in for the output of a real simulation and for a measured spectrum. A real fit
reads the files of the simulation instead. The script `make_files.jl` in the same
directory makes them again.

- There is one file for each simulated source, for example `cables_K40.lh5`, with the
  histograms `pdf/mul1/icpc` and `pdf/mul1/bege` (the two detector types) and the
  number of simulated primaries. The 2νββ decay has one file for each detector type.
- `data.lh5` has the measured histograms `hist/mul1/icpc` and `hist/mul1/bege`.

The sources with a γ line have a peak with a 1 keV sigma on a flat or linear
continuum below the peak. The 2νββ spectrum is a broad Gaussian (600 keV width,
centered at 1 MeV), and it makes about 70% of the counts. The ``^{39}``Ar spectrum is an
exponential with a decay constant of 300 keV and makes about 20% of the counts. The measured spectrum is
generated from known activities with Poisson fluctuations, so the result of the fit
can be compared with the truth at the end.

```@example tutorial
dir = "../../test/test_files/tutorial"  # relative to docs/src
readdir(dir)
```

## Settings

```@example tutorial
using Refit
using StatsBase

livetime = 1.0e6  # s
e_low, e_high = 100, 4000  # keV, the range of the fit
binning = 20  # merge 20 bins of 1 keV into one
blind = (2015.0, 2065.0)  # keV, a region around the 2039 keV Q-value that is not fitted

# expected counts = norm * activity * spectrum, see the note on units below
norm = 1e-3 * livetime

# the histograms of the measured spectrum and of each simulation are summed over these groups
spec = ["mul1/icpc", "mul1/bege"]
nothing # hide
```

!!! note "Units"
    In the standard case, where the PDFs are normalized, the fit parameters represent the
    activity of the sources. The unit of the activity is the inverse of the unit of the
    livetime, or Bq if the livetime is in s. In many cases it is more convenient to use a
    different unit. For example, scaling `norm` by 10⁻³ changes the unit of the activity
    to mBq, as done here.

## Reading the data

`read_hist` reads one histogram. The data of the two detector types are added with
`add_histograms`.

```@example tutorial
data = add_histograms([read_hist("hist/$s", joinpath(dir, "data.lh5")) for s in spec])
```

## Reading the simulated spectra

Each source has its own file. `read_pdf` sums the groups and divides by the number of
simulated primaries, so that the spectrum is in counts per primary. The name of each
spectrum is a `Symbol`, the activity of the source gets the same name.

```@example tutorial
groups = ["pdf/$s" for s in spec]

components = ["cables", "fibers"]
isotopes = ["K40", "Pb212_to_Pb208"]

mc_dict = Dict(
    Symbol(comp, "_", iso) => read_pdf(joinpath(dir, "$(comp)_$iso.lh5"); groups) for
    comp in components for iso in isotopes
)
keys(mc_dict)
```

The 2νββ decay is simulated once for each detector type. `read_pdf` accepts a
dictionary from file name to weight and returns the weighted average of the spectra.
The weights are normalised. Here the weight is the mass of each detector type in kg.

```@example tutorial
mass = Dict("icpc" => 80.0, "bege" => 20.0)
bulk = Dict(joinpath(dir, "hpge_$(t)_bulk_2vbb.lh5") => w for (t, w) in mass)

mc_dict[:two_vbb_hpge_bulk] = read_pdf(bulk; groups)
nothing # hide
```

The ``^{39}``Ar decay has one file.

```@example tutorial
mc_dict[:Ar39_liquid_argon] = read_pdf(joinpath(dir, "liquid_argon_Ar39.lh5"); groups)
nothing # hide
```

## Tying two activities

The activity of a component can also be a function of the fit parameters (see
[Activity functions](@ref)). As an example, assume that
the cables and the fibers are made of the same material, with the same ``^{40}``K
content per kg, and that the fibers have half the mass of the cables. Then the ``^{40}``K activity of the fibers is
half that of the cables, and needs no fit parameter of its own:

```@example tutorial
mass_ratio = 0.5  # fibers / cables

components = merge(
    mc_dict,
    Dict(:fibers_K40 => (pdf = mc_dict[:fibers_K40], activity = p -> mass_ratio * p.cables_K40)),
)
nothing # hide
```

The global `norm` multiplies this component like the others, so the function returns an
activity in mBq.

## Alpha events

Alpha events have no simulation. As a simple model they are a straight line in energy,
an analytic component (see [Analytic functions](@ref)). Its two parameters are the
counts per keV at the two ends of the fit range, so the line stays positive for any
values in their ranges. `norm` does not apply to it:

```@example tutorial
# from `alpha_low` counts per keV at 100 keV to `alpha_high` at 4000 keV
components[:alpha] = (
    analytic = (E, p) -> p.alpha_low + (p.alpha_high - p.alpha_low) * (E - 100) / 3900,
)
nothing # hide
```

The ends of the range are written as numbers: `e_low` and `e_high` are ordinary global
variables, which would make every call of the model much slower inside an analytic
function.

## Forward model and likelihood

The [forward model](forward_model.md) maps the activities to the expected counts in
each energy bin. The [likelihood](likelihood.md) compares these with the measured
counts, which are taken from the data in the same bins.

```@example tutorial
model, bin_edges = build_forward_model(components; e_low, e_high, blind, binning, norm)
counts = get_counts(data, bin_edges)
likelihood = build_likelihood(counts, model)
nothing # hide
```

## Ranges for the activities

The fit needs a range for each activity. `get_priors` estimates it from the data
(see [Priors](priors.md)). The reference ranges give the energy window of a
peak of each isotope in keV. The 1461 keV line of ``^{40}``K and the 2615 keV line of
``^{208}``Tl limit the activities of these sources. `get_priors` skips `fibers_K40`,
which has no parameter of its own, and the analytic `alpha` component, whose ranges
are given by hand.

```@example tutorial
using IntervalSets

ref_ranges = Dict(:K40 => (1456, 1466), :Pb212_to_Pb208 => (2610, 2620))

ranges = get_priors(components, data; ref_ranges, e_low, e_high, blind, norm)
ranges[:alpha_low] = 0.0 .. 10.0  # counts per keV
ranges[:alpha_high] = 0.0 .. 10.0
ranges
```

Any other range can be changed by hand in the same way, where the guess is poor.

`build_prior_dist` makes the prior: a uniform distribution on the range of each activity.
`rand(prior)` draws one random value for each activity, a named tuple that goes
directly into the model.

```@example tutorial
prior = build_prior_dist(ranges)

rand(prior)
```

## Sampling the posterior

[BAT.jl](https://bat.github.io/) samples the posterior with a Markov chain Monte
Carlo.

```@example tutorial
using BAT

posterior = PosteriorMeasure(likelihood, prior)

samples = bat_sample(
    posterior,
    TransformedMCMC(
        proposal = RandomWalk(),
        strict = false,
        burnin = MCMCMultiCycleBurnin(nsteps_per_cycle = 50_000, max_ncycles = 50),
        nsteps = 5 * 10^5,
        nchains = 4,
    ),
).result
nothing # hide
```

## Result

The activities used to make the fake data are below. The ``^{40}``K activity of the fibers
is half that of the cables, as assumed above.

```@example tutorial
true_activities = (  # mBq
    cables_K40 = 200.0,
    fibers_K40 = 100.0,
    cables_Pb212_to_Pb208 = 50.0,
    fibers_Pb212_to_Pb208 = 30.0,
    two_vbb_hpge_bulk = 3500.0,
    Ar39_liquid_argon = 1000.0,
)
true_alpha = 2.0  # counts per keV, the same at all energies
nothing # hide
```

The plots are made with CairoMakie, here in the LEGEND style of
[LegendMakie.jl](https://legend-exp.github.io/LegendMakie.jl/stable/), which changes the
style of all CairoMakie plots once it is loaded (it also needs `Format`). The best fit is
the global mode of the posterior. `plot_posterior!` draws the marginal posterior of one
parameter. The bins inside the smallest 68% interval are filled (`bandcolor` changes the color of the
band), the dashed line is the best fit and the red line is the true activity that was used
to make the fake data. All the panels have the same labels, so one legend is enough.

```@example tutorial
using LegendMakie
import Format
import CairoMakie

best = mode(samples)
truth = merge(true_activities, (alpha_low = true_alpha, alpha_high = true_alpha))

fig = CairoMakie.Figure(size = (1300, 700))
for (i, name) in enumerate(keys(best))
    unit = startswith(string(name), "alpha") ? "counts per keV" : "mBq"
    ax = CairoMakie.Axis(
        fig[fld1(i, 4), mod1(i, 4)];
        title = string(name),
        xlabel = "Activity [$unit]",
        ylabel = "Probability density",
        xticks = CairoMakie.LinearTicks(3),
    )
    plot_posterior!(
        ax,
        samples,
        name;
        best = best[name],
        truth = truth[name],
        bandcolor = :lightblue,
    )
end

# the legend of the first panel, in the empty place of the grid
CairoMakie.Legend(fig[2, 4], fig.content[1]; framevisible = false)
fig
```

The uncertainty of each activity is the smallest 68% credible interval, the smallest
interval that contains 68% of its posterior density, which is the filled region above.
BAT.jl computes it for all the activities at once, and `mode = :connected` gives a single
interval for each one. Activities differ by orders of magnitude, so the plot shows the
best fit divided by the true value. A good fit is consistent with 1.

```@example tutorial
intervals = map(first, BAT.smallest_credible_intervals(samples; p = 0.68, mode = :connected))

names = filter(in(keys(best)), collect(keys(truth)))  # the fit parameters
x = 1:length(names)
fit = [best[n] / truth[n] for n in names]
lower = [intervals[n].left / truth[n] for n in names]
upper = [intervals[n].right / truth[n] for n in names]

fig = CairoMakie.Figure(size = (900, 450))
ax = CairoMakie.Axis(
    fig[1, 1];
    xticks = (x, string.(names)),
    xticklabelrotation = pi / 6,
    ylabel = "Best fit / true activity",
)
CairoMakie.hlines!(ax, 1; color = :red, label = "True activity")
CairoMakie.errorbars!(ax, x, fit, max.(fit .- lower, 0), max.(upper .- fit, 0); color = :black)
CairoMakie.scatter!(ax, x, fit; color = :black, label = "Best fit and 68% interval")
CairoMakie.axislegend(ax; position = :rb)
fig
```

Finally, the prediction of the best fit is compared with the measured counts. The
sources are summed into a few groups and plotted in the same way. See
[Plotting](plotting.md) for the functions that are used and how the plot can be customised.

```@example tutorial
data_hist = counts_to_histogram(counts, bin_edges)
hists = get_best_fit_histograms(best, model, bin_edges)

# sum the sources into a few lines
grouped = group_histograms(
    hists,
    [
        "⁴⁰K" => r"K40",
        "²¹²Pb → ²⁰⁸Pb" => r"Pb212",
        "2νββ" => r"two_vbb",
        "³⁹Ar" => r"Ar39",
        "other" => :rest,
    ],
)

fig, ax, ax_res = fit_axes(; ylabel = "Counts per $binning keV")

# the data, the groups of sources (with the colors of the palette) and the best fit
plot_hist!(ax, data_hist; filled = true, alpha = 0.6, color = :lightblue, label = "Data")
for (label, h) in grouped
    label == :total || plot_hist!(ax, h; label = string(label))
end
plot_hist!(ax, grouped[:total]; color = :black, linewidth = 4, label = "Best fit")

# the blinded region (the empty bins) and the residuals
mark_blinded!(ax, grouped[:total])
mark_blinded!(ax_res, grouped[:total])
plot_residuals!(ax_res, grouped[:total], data_hist)

CairoMakie.xlims!(ax, e_low, e_high)
CairoMakie.ylims!(ax, 1, 1e5)
CairoMakie.ylims!(ax_res, -6, 6)
CairoMakie.axislegend(ax; nbanks = 3)
CairoMakie.axislegend(ax_res; nbanks = 3)
fig
```

## Saving the results

The results are saved in formats for long term preservation (see
[Saving the results](saving.md)). The posterior samples are written by BAT.jl to an
HDF5 file, which needs the package `HDF5`. The best fit histograms, of the total and of
each source, and the measured spectrum go into an LH5 file. The values of the activities
with the smallest 68% interval go into a YAML file.

```@example tutorial
using HDF5

outdir = mktempdir()

bat_write(joinpath(outdir, "samples.h5"), samples)
write_best_fit(joinpath(outdir, "best_fit.lh5"), hists; data = data_hist)
write_parameter_values(joinpath(outdir, "parameters.yaml"), best, intervals; level = 0.68)

readdir(outdir)
```

The YAML file is plain text:

```@example tutorial
print(read(joinpath(outdir, "parameters.yaml"), String))
```

The saved samples and histograms read back as they were written:

```@example tutorial
samples_saved = bat_read(joinpath(outdir, "samples.h5")).result
fit_saved = read_hist("best_fit/total", joinpath(outdir, "best_fit.lh5"))

(length(samples_saved) == length(samples), fit_saved.weights == hists[:total].weights)
```
