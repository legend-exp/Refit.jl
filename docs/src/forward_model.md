# Forward model

`build_forward_model` turns a set of simulated energy spectra into a function that
maps fit parameters to the predicted counts in each energy bin.

## Building the model

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
dictionary keys:

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

See the [API](api.md) for the full docstring.
