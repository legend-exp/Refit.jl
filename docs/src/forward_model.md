# Forward model

The forward model is a function mapping the parameters of interest onto the predicted
counts. Here the parameters are the activities of the simulated sources, for example
the decay rate of one isotope in one detector component. The prediction is the
expected number of counts in each energy bin of the measured spectrum.

`build_forward_model` turns a set of simulated energy spectra into such a function.
It maps the activities to the predicted counts in each energy bin.

## Building the model

This is a simplified example with two simulated spectra. A real fit has many more,
one for each source and each location in the setup. The spectrum names and the
file names are only examples.

```julia
using Refit, StatsBase

# simulated spectra (1D histograms, energy in keV) with identical binning
mc_dict = Dict(
    :Ar39_liquid_argon => read_hist("hist/spectrum", "Ar39_liquid_argon.lh5"),
    :two_vbb_hpge_bulk => read_hist("hist/spectrum", "2vbb_hpge_bulk.lh5"),
)

model, bin_edges = build_forward_model(mc_dict; e_low = 500, e_high = 3200)
```

The keys of `mc_dict` must be `Symbol`s. They name the fit parameters.
All histograms must have the same binning.

| Keyword | Default | Meaning |
|---|---|---|
| `e_low`, `e_high` | 500, 4000 | Energy range in keV. Only bins fully inside are kept |
| `blind` | `nothing` | `(min, max)` in keV. Bins overlapping this range are removed |
| `binning` | 1 | Number of consecutive bins merged before the range is applied |
| `norm` | 1 | Factor multiplying every spectrum, see below |

## Calling the model

The parameters are a named tuple with one activity per spectrum, named after the
dictionary keys. The activities here are made-up numbers, as in the example above
with two spectra:

```julia
params = (activities = (Ar39_liquid_argon = 1000.0, two_vbb_hpge_bulk = 0.5),)
pred = model(params)
```

`pred` has one entry per kept bin: the sum over all spectra of
`activity * norm * (counts in the simulated bin)`. Use `norm` for unit conversions
and the livetime, for example `norm = 1e-3 * livetime`. The order of the
activities does not matter. A missing activity raises an error.

## Bin edges

`bin_edges = (low = ..., high = ...)` holds the lower and upper edge in keV of each
kept bin, in the same order as `pred`. With `blind` the kept bins are not
contiguous, so two vectors are returned instead of one edges vector.

## Back to a histogram

`counts_to_histogram` turns the prediction into a `Histogram`, for example to plot
the model at the best-fit parameters. It is the inverse of `get_counts`
(see [Likelihood](likelihood.md)).

```julia
best_fit = counts_to_histogram(model(best_params), bin_edges)
```

The histogram covers the fit range, from `bin_edges.low[1]` to `bin_edges.high[end]`.
A blinded region inside the range becomes one bin with 0 counts.

See the [API](api.md) for the full docstring.
