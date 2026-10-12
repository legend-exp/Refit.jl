# Priors

`get_priors` guesses a range `0 .. a_max` for the activity of each simulated spectrum
from the measured spectrum.

!!! note
    These are not real priors. They come from the data and only limit the range of
    the fit to improve convergence. Use external knowledge, for example a screening
    measurement, for informative priors.

## The maximum activity

`a_max` is the largest activity for which the prediction does not exceed the data.
In a window between `low` and `high` keV it is

```math
a_\mathrm{max} = \frac{f \, N_\mathrm{data}}{\mathrm{norm} \, N_\mathrm{pdf}},
```

where ``N_\mathrm{data}`` and ``N_\mathrm{pdf}`` are the counts in the window of the
data and of the simulated spectrum, `norm` is the factor of the
[forward model](forward_model.md), and ``f`` is a safety factor. Only bins fully
inside the window and outside the blinded region are used.

Two windows are used and the smaller `a_max` is kept:
- the fit range `[e_low, e_high]`, with the safety factor `range_factor` (default 1.1),
- a reference range of the source, for example around a peak, with the safety factor
  `peak_factor` (default 1.5).

## Usage

```@setup demo
include("example_data.jl")
```

```@example demo
ref_ranges = Dict(
    :K40 => (1456, 1466),            # energy ranges in keV, here around a peak
    :Pb212_to_Pb208 => (2610, 2620),
)

priors = get_priors(
    mc_dict,
    data;
    ref_ranges,
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,  # as in build_forward_model
)
```

The result is a `Dict` with the same keys as `mc_dict`. Each value is an interval
`0 .. a_max`.

A spectrum uses the reference range of the longest key that its name starts with.
For example `:K40_hpge_bulk` uses `:K40`. Spectra without a range only use the fit
range. `e_low`, `e_high`, `blind` and `norm` are as in the forward model.

## Changing a range by hand

The result is a normal dictionary, so single entries can be replaced:

```@example demo
using IntervalSets

priors[:K40_hpge_bulk] = 0.0 .. 100.0
priors
```

Informative priors can use the guessed maximum:

```julia
using Distributions

a_max = priors[:Pb212_to_Pb208_hpge_bulk].right
priors[:Pb212_to_Pb208_hpge_bulk] = truncated(Normal(2.0, 0.5), 0, a_max)
```
