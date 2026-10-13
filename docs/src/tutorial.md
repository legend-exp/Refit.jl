# Tutorial: fitting a background spectrum

This tutorial fits a measured energy spectrum with a sum of simulated spectra, from
reading the files to sampling the posterior. It runs every time the documentation
is built.

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

The sources with a gamma line have a peak with a 1 keV width on a flat or linear
continuum below the peak. The 2νββ spectrum is a broad Gaussian (600 keV width,
centered at 1 MeV), and it makes about 90% of the counts. The measured spectrum is
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

# The activities are in mBq and the simulated spectra are counts per simulated
# primary, so the expected counts are 1e-3 * livetime * activity * spectrum
norm = 1e-3 * livetime

# the histograms of the measured spectrum and of each simulation are summed over these groups
spec = ["mul1/icpc", "mul1/bege"]
nothing # hide
```

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

Alpha events have no simulation. As a simple model they are a constant spectrum. Its
activity is the number of counts in each 1 keV bin, so the spectrum is
`1 / norm` in every bin.

```@example tutorial
mc_dict[:alpha] = Histogram(data.edges, fill(1 / norm, length(data.weights)))
nothing # hide
```

## Forward model and likelihood

The [forward model](forward_model.md) maps the activities to the expected counts in
each energy bin. The [likelihood](likelihood.md) compares these with the measured
counts, which are taken from the data in the same bins.

```@example tutorial
model, bin_edges = build_forward_model(mc_dict; e_low, e_high, blind, binning, norm)
counts = get_counts(data, bin_edges)
likelihood = build_likelihood(counts, model)
nothing # hide
```

## Ranges for the activities

The fit needs a range for each activity. `get_priors` estimates it from the data
(see [Priors](priors.md)). The reference ranges give the energy window of a
peak of each isotope in keV. The 1461 keV line of ``^{40}``K and the 2615 keV line of
``^{208}``Tl limit the activities of these sources.

```@example tutorial
ref_ranges = Dict(:K40 => (1456, 1466), :Pb212_to_Pb208 => (2610, 2620))

ranges = get_priors(mc_dict, data; ref_ranges, e_low, e_high, blind, norm)
```

A range can be changed by hand where the guess is poor, for example
`ranges[:two_vbb_hpge_bulk] = 0.0 .. 20.0` (with `using IntervalSets`).

`build_prior` makes the prior: a uniform distribution on the range of each activity.
`rand(prior)` draws one random value for each activity, a named tuple that goes
directly into the model.

```@example tutorial
prior = build_prior(ranges)

rand(prior)
```

## Sampling the posterior

[BAT.jl](https://bat.github.io/) samples the posterior with a Markov chain Monte
Carlo (a random walk through the space of the parameters that visits each value
with a probability given by the posterior). The first steps of each chain are
discarded as burn-in. This run is short, a real fit uses more steps.

```@example tutorial
using BAT

posterior = PosteriorMeasure(likelihood, prior)

samples = bat_sample(
    posterior,
    TransformedMCMC(
        proposal = RandomWalk(),
        strict = false,
        burnin = MCMCMultiCycleBurnin(nsteps_per_cycle = 10_000, max_ncycles = 30),
        nsteps = 10^5,
        nchains = 4,
    ),
).result
nothing # hide
```

## Result

The activities used to make the fake data are:

```@example tutorial
true_activities = (  # mBq
    cables_K40 = 200.0,
    fibers_K40 = 100.0,
    cables_Pb212_to_Pb208 = 50.0,
    fibers_Pb212_to_Pb208 = 30.0,
    two_vbb_hpge_bulk = 3500.0,
)
true_alpha = 2.0  # counts per 1 keV bin
nothing # hide
```

The marginal posterior of each activity is plotted with the plot recipes of
BAT.jl. The colored regions are the smallest intervals that contain 68%, 95% and
99.7% of the distribution, the dashed line is the mode and the red line is the true
activity used to make the measured spectrum.

```@example tutorial
using Plots

truth = merge(true_activities, (alpha = true_alpha,))
panels = map(collect(keys(truth))) do name
    unit = name == :alpha ? "counts per keV" : "mBq"
    plot(
        samples,
        name;
        mean = false,
        std = false,
        globalmode = true,
        nbins = 50,
        title = string(name),
        titlefontsize = 8,
        xlabel = "activity ($unit)",
        ylabel = "probability density (1/$unit)",
    )
    vline!([truth[name]]; color = :red, linewidth = 2, label = "")
end
plot(panels...; layout = (2, 3), size = (900, 650), legend = false)
```

A pair of activities shows their correlation. The K40 activities of the cables and of
the fibers are anticorrelated, because the data constrain mostly their sum:

```@example tutorial
plot(
    samples,
    (:cables_K40, :fibers_K40);
    nbins = 50,
    globalmode = true,
    legend = :bottomleft,
    xlabel = "cables_K40 activity (mBq)",
    ylabel = "fibers_K40 activity (mBq)",
)
```

The best fit is the global mode. The uncertainty of each activity is the central 68%
interval of its distribution. Both are compared with the activities used to make the
measured spectrum:

```@example tutorial
best = mode(samples)
weights = Weights(samples.weight)

println(rpad("source", 25), "best fit [68% interval]    true (mBq)")
for name in keys(true_activities)
    low, high = quantile([v[name] for v in samples.v], weights, [0.16, 0.84])
    println(
        rpad(name, 25),
        round(best[name], digits = 1),
        " [",
        round(low, digits = 1),
        ", ",
        round(high, digits = 1),
        "]    ",
        true_activities[name],
    )
end
low, high = quantile([v.alpha for v in samples.v], weights, [0.16, 0.84])
println(
    rpad("alpha (counts per keV)", 25),
    round(best.alpha, digits = 2),
    " [",
    round(low, digits = 2),
    ", ",
    round(high, digits = 2),
    "]    ",
    true_alpha,
)
```

Finally, the prediction of the best fit is compared with the measured counts, both as
histograms in the bins of the fit. `counts_to_histogram` turns the counts back into a
`Histogram`. The bins that overlap the blinded region are not fitted and are merged into one empty bin, shown as a grey band.

```@example tutorial
data_hist = counts_to_histogram(counts, bin_edges)
fit_hist = counts_to_histogram(model(best), bin_edges)

plot(
    fit_hist;
    seriestype = :stepbins,
    label = "best fit",
    color = :red,
    linewidth = 3,
    xlims = (e_low, e_high),
    xlabel = "energy (keV)",
    ylabel = "counts per $binning keV",
)
plot!(data_hist; seriestype = :stepbins, label = "data", color = :black, linewidth = 1)

# grey band over the empty bin, which is wider than the blinded region because bins
# are only fitted if they are fully outside of it
empty = findfirst(==(0), fit_hist.weights)
vspan!(fit_hist.edges[1][empty:(empty+1)]; color = :gray, alpha = 0.3, label = "blinded")
```
