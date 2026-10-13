# Model components

A fit model is a dictionary with one entry per component, each a simulated energy
spectrum:

```julia
Dict(
    :Bi212_cables => h_Bi212,                                        # a fit parameter
    :Tl208_cables => (pdf = h_Tl208, activity = p -> 0.3594 * p.Bi212_cables),
)
```

The key is a `Symbol` and names the component. The value gives its spectrum and says
where its activity comes from. This page lists the forms a value can take and the rules
that follow from them.

## The two forms

| Value | Activity ``A_k`` of the component |
|---|---|
| `h` | `params[name]`, the fit parameter with the name of the component |
| `(pdf = h, activity = f)` | `f(params)`, a function of the fit parameters |

`h` is a `Histogram` of a simulated 1D energy spectrum in keV and `params` the named
tuple of **all** fit parameters. `(pdf = h,)` is the same as `h`. Any other key, such
as a misspelled `activty`, is an error.

The predicted counts in bin ``i`` are the sum over the components:

```math
\mu_i = \sum_k A_k \times N \times h_{k,i}
```

with ``h_{k,i}`` the content of bin ``i`` of spectrum ``k`` and ``N`` the global `norm`
of [`build_forward_model`](@ref), usually the unit conversion times the livetime. It
applies to every component, so an `activity` function returns an activity in the same
units as a fit parameter, not counts.

## Parameters

A component of the first form takes its parameter from its own name: the entry
`:Bi212_cables` needs `params.Bi212_cables`.

An `activity` function takes whatever it reads out of `params`, under any name, and
several components can read the same parameter. Above, the activity of Tl-208 is tied
to that of Bi-212 by the branching ratio of 35.94 %, so `:Tl208_cables` adds no
parameter of its own. A function can also replace the activity by the quantity of
interest, for example a half-life, see [Activity functions](forward_model.md#Activity-functions).

A parameter missing from `params` raises an error at the first call of the model,
naming the field it did not find. Parameters not used by any component are ignored.

## Priors

[`get_priors`](@ref) estimates a range only for the components whose activity is a fit
parameter. It skips the others with a message, and the parameters of their functions
need priors by hand:

```julia
priors = get_priors(components, data; e_low = 500, e_high = 3200, norm = 1e-3 * livetime)
priors[:T_half] = 1.0 .. 5.0
```

## Each component separately

`model(params, Val(:by_component))` returns a named tuple with the predicted counts of
each component on its own, under the keys of the dictionary and in sorted order,
whichever form the component has. [`get_best_fit_histograms`](@ref) uses this to draw
the fit.
