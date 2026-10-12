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
- the fit range `[e_low, e_high]`, with `range_factor` (default 1.1),
- a window around a reference peak of the source, with `peak_factor` (default 1.5).

## Usage

```julia
ref_regions = Dict(
    :K40 => 1461,                    # a peak at 1461 keV, the window is 1461 ± 5 keV
    :Pb212_to_Pb208 => 2615,
    :Ra228_to_Ac228 => (1580, 1596), # a window in keV
)

priors = get_priors(
    mc_dict,
    data;
    ref_regions,
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,  # as in build_forward_model
)
```

The result is a `Dict` with the same keys as `mc_dict`. Each value is an interval
`0 .. a_max`.

A spectrum uses the reference region of the longest key that its name starts with.
For example `:K40_hpge_bulk` uses `:K40`. Spectra without a region only use the fit
range. `peak_width` (default 5 keV), `peak_factor`, `range_factor`, `e_low`,
`e_high`, `blind` and `norm` can all be changed.

## Changing a range by hand

The result is a normal dictionary, so single entries can be replaced. Informative
priors can use the guessed maximum:

```julia
priors[:K40_hpge_bulk] = 0.0 .. 3.0

using Distributions
a_max = priors[:Pb212_to_Pb208_hpge_bulk].right
priors[:Pb212_to_Pb208_hpge_bulk] = truncated(Normal(2.0, 0.5), 0, a_max)
```
