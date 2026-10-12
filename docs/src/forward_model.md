# Forward model

The forward model is a function mapping the parameters of interest onto the predicted
counts. Here the parameters are the activities of the simulated sources, for example
the decay rate of one isotope in one detector component. The prediction is the
expected number of counts in each energy bin of the measured spectrum.

`build_forward_model` turns a set of simulated energy spectra into such a function.
It maps the activities to the predicted counts in each energy bin.

## Building the model

The examples in this manual use made-up spectra with three simulated sources, defined
in `example_data.jl`. A real fit has many more, one for each source and each location
in the setup. The simulated spectra are a dictionary of `Histogram`s with the energy
in keV, and all have the same binning:

```@setup demo
include("example_data.jl")
```

```@example demo
keys(mc_dict)
```

```@example demo
model, bin_edges = build_forward_model(mc_dict; e_low = 500, e_high = 3200)
```

`model` is the function mapping parameters to predicted counts, `bin_edges` describes
its bins (see below).

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
dictionary keys. The activities here are in mBq, so `norm` must convert them to
counts:

```@example demo
model, bin_edges = build_forward_model(
    mc_dict;
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,  # mBq to Bq, times the livetime in seconds
)

params = (activities = (
    K40_hpge_bulk = 40.0,
    Pb212_to_Pb208_hpge_bulk = 10.0,
    two_vbb_hpge_bulk = 100.0,
),)

pred = model(params)
```

`pred` is a flat vector with one entry per kept bin. Entry `i` is the sum over all
spectra of `activity * norm * (counts in the simulated bin)`. The order of the
activities does not matter. A missing activity raises an error.

## Bin edges

`bin_edges` holds the lower and upper edge in keV of each kept bin, in the same order
as `pred`:

```@example demo
bin_edges.low[1:3], bin_edges.high[1:3]
```

```@example demo
length(pred) == length(bin_edges.low)
```

With `blind` the kept bins are not contiguous, so two vectors are returned instead of
one edges vector.

## Back to a histogram

`counts_to_histogram` turns the prediction into a `Histogram`, for example to plot
the model at the best-fit parameters. It is the inverse of `get_counts`
(see [Likelihood](likelihood.md)).

```@example demo
best_fit = counts_to_histogram(pred, bin_edges)
best_fit.edges[1][[1, end]], length(best_fit.weights)
```

The histogram covers the fit range, from `bin_edges.low[1]` to `bin_edges.high[end]`.
A blinded region inside the range becomes one bin with 0 counts.

See the [API](api.md) for the full docstring.
