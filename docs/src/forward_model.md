# Forward model

The forward model is a function mapping the parameters of interest onto the predicted
observables (in our case the counts in each bin).

## Definition
In `Refit.jl` a linear template fit model is implemented,
the predicted counts in energy bin ``i`` are a sum over the simulated spectra:

```math
\mu_i(\vec{a}) = \sum_k a_k \times h_{k,i} \times N
```

Here ``h_{k,i}`` is the content of bin ``i`` in simulated spectrum ``k``, ``N`` is the
factor `norm`, and ``a_k`` is the parameter of spectrum ``k``. Each simulated spectrum
is scaled by its parameter and by `norm`, and the scaled spectra are added.

If `norm` holds the unit conversions and the livetime, ``a_k`` is the activity of
source ``k``. Usually the parameters are activities, and the model is written for
that case. Any other choice of `norm` changes the meaning of ``a_k`` in the same way.

`build_forward_model` turns a set of simulated energy spectra into such a function.
This maps the activities to the predicted counts in each energy bin.

The activity of a component can also be a function of the fit parameters, see
[Activity functions](#Activity-functions) below.

## Building the model

The examples in this manual use simplified spectra with three simulated sources, defined
in `example_data.jl`. A real fit has many more sources.
The simulated spectra are a dictionary of `Histogram`s with the energy
in keV, and all have the same binning:

```@setup demo
include("example_data.jl")
```

```@example demo
# the names of the simulated spectra, these name the parameters
keys(mc_dict)
```

```@example demo
# keep the bins between 500 and 3200 keV
model, bin_edges = build_forward_model(
    mc_dict;
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,  # mBq to Bq, times the livetime in seconds
)
nothing # hide
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

The parameters are a flat named tuple with one activity per spectrum, named after the
dictionary keys. Here the activities are in mBq, which `norm` converts to counts:

```@example demo
params = (
    K40_hpge_bulk = 40.0,
    Pb212_to_Pb208_hpge_bulk = 10.0,
    two_vbb_hpge_bulk = 100.0,
)

# predicted counts in each bin
pred = model(params)
```

`pred` is a flat vector with one entry per kept bin, here
`length(pred)` = 270. Entry `i` is the sum over all
spectra of `activity * norm * (counts in the simulated bin)`. The order of the
activities does not matter. A missing activity raises an error.

## Activity functions

The activity of a component does not have to be a fit parameter of its own. It can be a
function of the fit parameters, given with the spectrum as `(pdf = histogram, activity
= f)`. For example, to fit the half-life of the ``2\nu\beta\beta`` decay instead of its
activity:

```@example demo
n_ge76 = 1.0e28  # made-up number of Ge-76 nuclei
year = 3.156e7   # s

components = merge(
    mc_dict,
    Dict(
        :two_vbb_hpge_bulk => (
            pdf = mc_dict[:two_vbb_hpge_bulk],
            # activity in mBq, with the half-life in units of 1e21 yr
            activity = p -> 1e3 * log(2) * n_ge76 / (p.T_half * 1e21 * year),
        ),
    ),
)
model_t, _ = build_forward_model(components; e_low = 500, e_high = 3200, norm = 1e-3 * livetime)
pred_t = model_t((K40_hpge_bulk = 40.0, Pb212_to_Pb208_hpge_bulk = 10.0, T_half = 2.2))

# the activity in mBq for a half-life of 2.2e21 yr
components[:two_vbb_hpge_bulk].activity((T_half = 2.2,))
```

The global `norm` applies to this component as to the others, so the function returns
an activity in mBq. [Model components](components.md) lists the forms a component can
take and the rules for its parameters and priors.

`model(params, Val(:by_component))` returns the predicted counts of each component
separately, as a named tuple. [`get_best_fit_histograms`](@ref) uses this.

## Bin edges

`bin_edges` holds the lower and upper edge in keV of each kept bin, in the same order
as `pred`:

```@example demo
# lower, upper edge of each bin in keV
bin_edges.low, bin_edges.high

```

With `blind` the kept bins are not contiguous, so two vectors are returned instead of
one edges vector.

## Back to a histogram

`counts_to_histogram` turns the prediction into a `Histogram`, for example to plot
the model at the best-fit parameters. It is the inverse of `get_counts`
(see [Likelihood](likelihood.md)).

```@example demo
# the prediction as a Histogram, for plotting
best_fit = counts_to_histogram(pred, bin_edges)
```

The histogram covers the fit range, from `bin_edges.low[1]` to `bin_edges.high[end]`.
A blinded region inside the range becomes one bin with 0 counts.

See the [API](api.md) for the full docstring.
