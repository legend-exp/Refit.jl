# Forward model

The forward model is a function mapping the parameters of interest onto the predicted
observables, in our case the counts in each bin. This defines the model, is
used in the likelihood, plotting etc.

## Definition

In `Refit.jl` the model is a sum of components. The predicted counts in energy bin ``i``
are

```math
\mu_i(\vec{a}) = \sum_k \mu_{k,i}(\vec{a}),
```

with ``\vec{a}`` the fit parameters. A component is one of two kinds:

- a simulated spectrum scaled by an activity, a template:

  ```math
  \mu_{k,i}(\vec{a}) = A_k(\vec{a}) \times N \times g_k(E_i, \vec{a}) \times h_{k,i}
  ```

- an analytic spectrum, with no simulation:

  ```math
  \mu_{k,i}(\vec{a}) = g_k(E_i, \vec{a}) \times w_i
  ```

Here ``h_{k,i}`` is the content of bin ``i`` in the simulated spectrum of component ``k``,
``N`` is the factor `norm`, and ``E_i`` and ``w_i`` are the energy at the centre and the
width of bin ``i`` in keV.
- ``A_k`` is the activity of the component. By default it is a fit parameter of its own,
  ``A_k = a_k``. It can also be a function of the fit parameters, ``A_k = f_k(\vec{a})``,
  for example to tie the activities of two components together or to fit a half-life
  instead of an activity (see [Activity functions](@ref) below).
- ``g_k`` is a function of the energy and of the fit parameters (see
  [Analytic functions](@ref) below). With a simulated spectrum it multiplies the
  spectrum bin by bin, for example an efficiency curve, and ``g_k = 1`` by default. For
  an analytic spectrum it is the expected number of counts per keV.

``N`` multiplies every component with a simulated spectrum, also those whose activity is
a function, but not an analytic component with no simulation. If `norm` holds the unit
conversions and the livetime, ``A_k`` is the activity of source ``k``, and an activity
function returns an activity in the same units as a fit parameter. Usually the
parameters are activities, and the model is written for that case. Any other choice of
`norm` changes the meaning of ``A_k`` in the same way.

`build_forward_model` turns a set of components into such a function. This maps the
parameters to the predicted counts in each energy bin. It is one way to make a forward
model: any function of the parameters that returns the predicted counts in the bins of
the measured counts can be used in the [likelihood](likelihood.md), also one written by
hand.

## Model components

The model is given as a dictionary with one entry per component. The key is a `Symbol`
and names the component. The value says how its predicted counts are built:

| Value | Predicted counts in bin ``i`` |
|---|---|
| `h` | `params[name] * norm * h[i]` |
| `(pdf = h, activity = f)` | `f(params) * norm * h[i]` |
| `(pdf = h, analytic = g)` | `params[name] * norm * g(E_i, params) * h[i]` |
| `(pdf = h, activity = f, analytic = g)` | `f(params) * norm * g(E_i, params) * h[i]` |
| `(analytic = g,)` | `g(E_i, params) * w_i` |

`h` is a `Histogram` of a simulated 1D energy spectrum in keV, `h[i]` its content in bin
``i``, and `params` the named tuple of **all** fit parameters. All histograms must have
the same binning. `(pdf = h,)` is the same as `h`. An analytic component with no `pdf`
has no activity. Any other key, such as a misspelled `activty`, is an error.

`h` can also be a [`GeneralisedHistogram`](@ref), a simulated spectrum that depends on
fit parameters, such as the thickness of a dead layer (see
[Accounting for systematic uncertainties](generalised_hist.md)). Then `h[i]` is
interpolated at the values of these parameters in `params` at every evaluation.

A component with a `pdf` and no `activity` function takes its parameter from its own
name: the entry `:K40_hpge_bulk` needs `params.K40_hpge_bulk`. An activity or analytic
function takes whatever it reads out of `params`, under any name, and several
components can read the same parameter. A parameter missing from `params` raises an
error at the first call of the model, naming the field it did not find. Parameters not
used by any component are ignored.

## Building the model

The examples in this manual use simplified spectra with three simulated sources, defined
in `example_data.jl`.
The simulated spectra are a dictionary of `Histogram`s with the energy
in keV, and all have the same binning:

```@setup demo
include("example_data.jl")
```

```@example demo
# the names of the components, here each also names a fit parameter
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

| Keyword | Default | Meaning |
|---|---|---|
| `e_low`, `e_high` | 500, 4000 | Energy range in keV. Only bins fully inside are kept |
| `binning` | 1 | Number of consecutive bins merged before the range is applied |
| `norm` | 1 | Factor ``N`` multiplying every component with a simulated spectrum |
| `edges` | `nothing` | Bin edges in keV, needed only for a model with no histogram |
| `points_per_bin` | 1 | Number of points in each bin at which an analytic function is evaluated |

`binning` can also give bins of any width and leave regions out, for example a blinded
one, see [Binning](binning.md).

## Calling the model

Here every component is of the first form, so the parameters are a flat named tuple with
one activity per spectrum, named after the dictionary keys. The activities are in mBq,
which `norm` converts to counts:

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

The activity of a component does not have to be a fit parameter of its own. For example,
the activity of ``^{208}``Tl is tied to that of ``^{212}``Bi by the branching ratio of
35.94 %, so this component adds no parameter of its own:

```julia
:Tl208_cables => (pdf = h_Tl208, activity = p -> 0.3594 * p.Bi212_cables)
```

A function can also replace the activity by the quantity of interest. To fit the
half-life of the ``2\nu\beta\beta`` decay instead of its activity:

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

The parameters of this model are `K40_hpge_bulk`, `Pb212_to_Pb208_hpge_bulk` and
`T_half`. `norm` multiplies the 2νββ component like the others, so the function returns
an activity in mBq.

[`get_priors`](@ref) only estimates a range for activities that are fit parameters. It
skips the components with an activity function, and the parameters of their functions
need a prior by hand, see [Priors](priors.md).

## Analytic functions

An analytic function `g(E, params)` depends on the energy `E` in keV and on the fit
parameters. On its own, with no simulated spectrum, it is a spectrum: the expected
number of counts per keV, which is multiplied by the width of each bin. `norm` does not
apply to it. For example a flat background of 0.5 counts per keV:

```@example demo
components_g = merge(mc_dict, Dict(:flat => (analytic = (E, p) -> p.b,)))
model_g, _ = build_forward_model(components_g; e_low = 500, e_high = 3200, norm = 1e-3 * livetime)

# 0.5 counts per keV in bins of 10 keV
model_g(merge(params, (b = 0.5,)), Val(:by_component)).flat[1:3]
```

Together with a simulated spectrum, the function multiplies it bin by bin, on top of the
activity and `norm`. For example the spectrum of a beta decay
``S(E) = F(E) \, (Q - E)^2``, with ``F(E)`` from a simulation and the endpoint ``Q`` a fit
parameter, or an efficiency curve:

```julia
:beta => (pdf = h_F, analytic = (E, p) -> max(p.Q - E, 0)^2)
:K40_hpge_bulk => (pdf = h_K40, analytic = (E, p) -> 1 - exp(-E / p.E0))
```

A model can be made of analytic components only. There is then no histogram to take the
binning from, so the bin edges in keV are given with `edges`, see the tutorial
[Fitting a peak](tutorial_peak.md).

By default an analytic function is evaluated at the centre of each bin only. A feature
narrower than about one bin, such as a peak, then comes out wrong. With
`points_per_bin = n`, `build_forward_model` splits each bin into `n` equal parts,
evaluates the function at the centre of each and uses the average. This takes `n` times
as many calls of the function, and the error of the integral over the bin falls as
``1/n^2``: `points_per_bin = 5` makes it 25 times smaller.

An analytic function with no simulated spectrum must not be negative for any value of
the parameters inside their prior ranges: a negative prediction has no Poisson
probability, the log-likelihood is `-Inf` and these parameters are rejected, which
wastes part of the prior range. Write it so that it cannot be negative, for example a
straight line by its values at the two ends of the fit range.

!!! warning "Constants in an analytic function"
    An analytic function is called in every bin, on every call of the model. A constant
    used inside it, such as a livetime, must be a `const` or a global with a fixed type,
    for example `livetime::Float64 = 1e6`. An ordinary global variable makes every call
    of the model about ten times slower, because Julia cannot know its type when it
    compiles the function. `build_forward_model` warns when an analytic function reads
    one.

Like for activity functions, `get_priors` skips the analytic components, and their
parameters need a prior by hand.

!!! note "Each component separately"
    `model(params, Val(:by_component))` returns a named tuple with the predicted counts
    of each component on its own, under the keys of the dictionary and in sorted order,
    whichever form the component has. [`get_best_fit_histograms`](@ref) uses this to
    draw the fit.

## Bin edges

`bin_edges` holds the lower and upper edge in keV of each kept bin, in the same order
as `pred`:

```@example demo
# lower, upper edge of each bin in keV
bin_edges.low, bin_edges.high

```

With regions left out by the binning, such as a blinded one, the kept bins are not
contiguous, so two vectors are returned instead of one edges vector.

## Back to a histogram

`counts_to_histogram` turns the prediction into a `Histogram`, for example to plot
the model at the best-fit parameters. It is the inverse of `get_counts`
(see [Likelihood](likelihood.md)).

```@example demo
# the prediction as a Histogram, for plotting
best_fit = counts_to_histogram(pred, bin_edges)
```

The histogram covers the fit range, from `bin_edges.low[1]` to `bin_edges.high[end]`.
A region skipped by the binning, such as a blinded one, becomes one bin with 0 counts.

See the [API](api.md) for the full docstring.
