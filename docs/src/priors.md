# Priors

`get_priors` guesses a range `0 .. a_max` for the activity of each simulated spectrum
from the measured spectrum.

!!! note
    These are not real priors. They come from the data and only limit the range of
    the fit to improve convergence. Use external knowledge, for example a screening
    measurement, for informative priors (see
    [Priors from measured activities](#Priors-from-measured-activities)).

## The maximum activity

`a_max` is the largest activity for which the prediction does not exceed the data.
In a window between `low` and `high` keV it is

```math
a_\mathrm{max} = \frac{f \, N_\mathrm{data}}{\mathrm{norm} \, N_\mathrm{pdf}},
```

where ``N_\mathrm{data}`` and ``N_\mathrm{pdf}`` are the counts in the window of the
data and of the simulated spectrum, `norm` is the factor of the
[forward model](forward_model.md), and ``f`` is a safety factor. Only bins fully
inside the window are used, and, if `binning` has several ranges, only those inside one
of them: a blinded region left out of the fit is not looked at.

Two windows are used and the smaller `a_max` is kept:
- the fit range `[e_low, e_high]`, with the safety factor `range_factor` (default 1.1),
- a reference range of the source, for example around a peak, with the safety factor
  `peak_factor` (default 1.5).

!!! note "Generalised histograms"
    For a [generalised histogram](generalised_hist.md), whose spectrum depends on
    parameters, `a_max` is computed for the spectrum at every grid point, and the largest
    is used. The range then holds the activity for any values of the parameters in the
    grid, usually set by the grid point with the fewest simulated counts. When both
    windows are used, the widest range can lie slightly between grid points, so it is
    an estimate, like the other ranges.

## Usage

```@setup demo
include("example_data.jl")
```

```@example demo
# energy range in keV for each source, here around a peak
ref_ranges = Dict(
    :K40 => (1456, 1466),
    :Pb212_to_Pb208 => (2610, 2620),
)

# a range 0 .. a_max for the activity of each spectrum
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
range. `e_low`, `e_high`, `binning` and `norm` are as in the forward model.

## Changing a range by hand

The result is a dictionary, so single entries can be replaced:

```@example demo
using IntervalSets

# replace the range of one spectrum by hand
priors[:K40_hpge_bulk] = 0.0 .. 100.0
priors
```

A range can also come from another energy window with `get_max_activity`. It returns
the activity ``a_\mathrm{max}`` of the formula above without the safety factor, for any
window. Here the range of ``^{40}``K uses only the counts of its 1461 keV line, with a
safety factor of 1.5:

```@example demo
norm = 1e-3 * livetime  # as in build_forward_model
priors[:K40_hpge_bulk] = 0.0 .. 1.5 * get_max_activity(1450, 1470, data, mc_dict[:K40_hpge_bulk]; norm)
```

The prior does not have to be a range. It can be any 1D distribution from
[Distributions.jl](https://juliastats.org/Distributions.jl/stable/univariate/).
For example, a Gaussian with mean 50 and width 10, truncated to positive
activities between 0 and 1000.

```@example demo
using Distributions

priors[:K40_hpge_bulk] = truncated(Normal(50.0, 10.0), 0.0, 1000.0)
priors
```

## Priors from measured activities

A measured activity, for example from the screening of a material, gives an informative
prior. `activity_prior` makes a Gaussian cut at 0 from a value and its uncertainty, or a
half-normal distribution from an upper limit, whose quantile at the confidence level is
the limit. The activities have units (from Unitful.jl). The distribution is of the
activity as a number in `unit`, which must be the unit of the activities of the forward
model, here mBq:

```@example demo
using Unitful

priors[:K40_hpge_bulk] = activity_prior(; val = 42u"mBq", unc = 5u"mBq", unit = u"mBq")
activity_prior(; limit = 15u"mBq", cl = 0.9, unit = u"mBq")  # an upper limit at 90% CL
```

The half-normal distribution of an upper limit is highest at 0. `limit_prior` chooses
another one, with the same quantile at the confidence level: `:exponential`, also highest
at 0 but with more weight on low values, its median at 0.30 times the limit for 90% CL
instead of 0.41, or `:lognormal`, a log-normal distribution of width 1 in natural-log
units, which is 0 at 0 and has its median at 0.28 times the limit, or
[`lognormal_limit`](@ref) with another width, or any function of the limit and the
confidence level:

```@example demo
activity_prior(; limit = 15u"mBq", cl = 0.9, unit = u"mBq", limit_prior = :exponential)
activity_prior(; limit = 15u"mBq", cl = 0.9, unit = u"mBq", limit_prior = :lognormal)
activity_prior(; limit = 15u"mBq", cl = 0.9, unit = u"mBq", limit_prior = lognormal_limit(σ = 2))
activity_prior(; limit = 15u"mBq", cl = 0.9, unit = u"mBq", limit_prior = (L, cl) -> Uniform(0, L / cl))
```

The activities of the components of an experiment can be listed in a YAML file. Each
component has a quantity, its mass in kg, its area in m² or its number of pieces, and the
activity per unit of quantity of each simulation (`simid`) of its decays, in µBq per kg,
per m² or per piece. An activity is a measured value with its uncertainty, or an upper
limit with its confidence level:

```yaml
hpge_detectors:
  mass_in_kg: 20
  specific_activity_in_uBq_per_kg:
    K40_hpge_bulk: { val: 2100, unc: 250 }
    Pb212_to_Pb208_hpge_bulk: { limit: 750, cl: 0.9 }
```

```@setup demo
write(
    "activities.yaml",
    """
    hpge_detectors:
      mass_in_kg: 20
      specific_activity_in_uBq_per_kg:
        K40_hpge_bulk: { val: 2100, unc: 250 }
        Pb212_to_Pb208_hpge_bulk: { limit: 750, cl: 0.9 }
    """,
)
```

`read_activities` reads the file and gives the prior of the activity of each component
and simulation, the quantity times the activity per unit. `get_activity` takes the prior
of one simulation from it:

```@example demo
activities = read_activities("activities.yaml"; unit = u"mBq")
priors[:Pb212_to_Pb208_hpge_bulk] = get_activity(activities, :Pb212_to_Pb208_hpge_bulk)
```

A production cycle of legend-simflow has this file in its configuration, and
`read_activities(SimProd(path))` reads it. A simulation that stands for several
components, for example parts of different materials, has the sum of their activities:
`get_activity` then gives an error, and the prior is left to the user.

## Components with an activity or analytic function

`get_priors` can only estimate a range for an activity that is a fit parameter, of a
simulated spectrum that is not changed by an energy dependent function. It skips the
components with an `activity` or an `analytic` function, with a message, and the
parameters of their functions need a prior by hand. Here the half-life of the ``2\nu\beta\beta`` decay
replaces its activity (see [Activity functions](@ref)):

```@example demo
components = merge(
    mc_dict,
    Dict(
        :two_vbb_hpge_bulk => (
            pdf = mc_dict[:two_vbb_hpge_bulk],
            # activity in mBq, with the half-life in units of 1e21 yr
            activity = p -> 1e3 * log(2) * 1.0e28 / (p.T_half * 1e21 * 3.156e7),
        ),
    ),
)
priors_t = get_priors(
    components,
    data;
    ref_ranges,
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,
)
priors_t[:T_half] = 1.0 .. 5.0  # in units of 1e21 yr
priors_t
```

## The prior distribution

`build_prior_dist` turns the entries into a prior. A range becomes a uniform distribution
on that range, a distribution is used as it is. `rand(prior)` draws one value for each activity as a named tuple, which
goes directly into the forward model.

```@example demo
prior = build_prior_dist(priors)

rand(prior)
```
