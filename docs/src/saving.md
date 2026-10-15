# Saving the results

The results of a fit are saved for long term preservation and further analysis. The
posterior samples, the best fit spectra as histograms in LH5 files, and the values of the
parameters with their uncertainty in a YAML file. The [tutorial](tutorial.md) shows all of them
for a complete fit.

```@setup demo
include("example_data.jl")

using Refit

# the fit of the forward model page, with the true activities as the best fit
norm = 1e-3 * livetime
# 2010 to 2070 keV is blinded: it is left out between the two ranges of the binning
binning = [500:10:2010, 2070:10:3200]
model, bin_edges = build_forward_model(mc_dict; e_low = 500, e_high = 3200, binning, norm)
best = true_activities
outdir = mktempdir()
```

## Posterior samples

The samples are saved by [BAT.jl](https://bat.github.io/), which writes them to an
HDF5 file with `bat_write` and reads them back with `bat_read`. The package `HDF5`
has to be loaded:

```julia
using BAT, HDF5

bat_write("samples.h5", samples)
samples = bat_read("samples.h5").result
```

The file has one dataset for each quantity, with one entry for each sample:

```text
v/<name>        value of each parameter (here the activities)
weight          number of times the sample was repeated by the sampler
logd            logarithm of the posterior density at the sample
info/chainid    chain, walker, step number and other bookkeeping of the sampler
```

## Best fit spectra

`get_best_fit_histograms` makes a histogram of the prediction of the model for the
`parameters`, usually the best fit. The result is a dictionary with the total in the
key `:total` and the contribution of each source, with the name of its activity. The
bins outside the fit range and in a region skipped by the binning, such as a blinded
one, are empty.

```@example demo
hists = get_best_fit_histograms(best, model, bin_edges)

keys(hists)
```

`write_best_fit` writes them to an LH5 file, with the total in `best_fit/total` and
each source in `best_fit/<name>`. With `data` the measured spectrum in the same bins is
saved too, in `data`. All can be read back with `read_hist`.

```@example demo
file = joinpath(outdir, "best_fit.lh5")

counts = get_counts(data, bin_edges)
write_best_fit(file, hists; data = counts_to_histogram(counts, bin_edges))

read_hist("best_fit/K40_hpge_bulk", file)
```

## Parameter values

`write_parameter_values` writes the values of the parameters, here the activities, to a
YAML file. For each parameter it saves the best fit and an interval, and the file
records the credible level of the intervals. The intervals are the smallest credible
intervals: the shortest intervals that contain 68% of the posterior. BAT.jl computes them, and `mode = :connected` gives one
interval for each parameter:

```julia
using Statistics

best = mode(samples)
intervals = map(first, BAT.smallest_credible_intervals(samples; p = 0.68, mode = :connected))
```

!!! note
    BAT.jl marks `smallest_credible_intervals` as internal, so it may change
    between versions.

Here the intervals are made up:

```@example demo
using IntervalSets

file = joinpath(outdir, "parameters.yaml")

intervals = map(a -> 0.9a .. 1.1a, best)
write_parameter_values(file, best, intervals; level = 0.68)

print(read(file, String))
```
