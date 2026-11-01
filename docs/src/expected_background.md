# Expected background from known activities

The [forward model](forward_model.md) predicts the counts in each energy bin for one set
of parameters, for example the activities of the sources. Passing many samples of the
parameters through it gives the distribution of the expected spectrum and of any quantity
derived from it, such as an event rate or the background index at the Q-value of the 0νββ
decay. The samples can come from a fit to a measured spectrum, or be drawn from the priors
before there is one.

This page draws them from priors on the activities of the components of an experiment,
from radioassay measurements, upper limits or decay rates, to predict its background.

## Priors of the activities

[`activity_prior`](@ref) turns a measured activity, with its uncertainty, into a Gaussian
cut at 0, and an upper limit into a half-normal distribution whose quantile at the
confidence level is the limit. The activities have units (from Unitful.jl), and the
distribution is of the activity in µBq by default:

```julia
using Unitful
activity_prior(; val = 8.1u"mBq", unc = 1.2u"mBq")  # a measurement
activity_prior(; limit = 160u"μBq", cl = 0.9)        # an upper limit at 90% confidence
```

[`read_activities`](@ref) reads these priors from a YAML file, which lists the components
of the experiment, their quantity (mass, area or number of pieces) and the activity per
unit of quantity of each simulation (`simid`) that stands for them. The prior of a
component is the quantity times the activity per unit. A production cycle of
legend-simflow has the file in its configuration directory:

```julia
prod = SimProd("prod/l1000dsg01/v1.1.0")
activities = read_activities(prod)  # Dict(component => Dict(simid => prior))
```

The result is a dictionary, so priors can be added or changed by hand, for example with an
activity from a half-life:

```julia
# 2νββ decay of the 7.3e27 nuclei of 76Ge in 1000 kg of germanium enriched to 92%, with a
# half-life of 2.02e21 yr, with an uncertainty of 2%
a = log(2) / 2.02e21u"yr" * 7.3e27
activities[:two_nu_bb] = Dict(:hpge_icpc_bulk_2vbb => activity_prior(; val = a, unc = 0.02 * a))
```

## Sources of the forward model

[`activity_components`](@ref) makes one source of the forward model for each component
and simid, named `<component>_<simid>`, with the spectrum of the simid and the prior of
the activity of the component. A simulation that stands for parts of different materials
gives one source for each part. With the production cycle, the spectra are read with the
event selection of [`read_pdf`](@ref):

```julia
components, priors = activity_components(prod, activities; groups = [:is_ss_aoe, :not_spms_coincident])
```

`norm` turns an activity in µBq into the number of decays in the live time. With the
live time as a quantity with units, it is the number µBq × live time:

```julia
livetime = 10u"yr"
norm = ustrip(NoUnits, u"μBq" * livetime)
model, bin_edges = build_forward_model(components; norm, e_low = 500, e_high = 3000, binning = 10)
```

The activities are drawn from their priors. The samples have the format of those of
BAT.jl, so the samples of a fit can be used in the same way below:

```julia
prior = build_prior_dist(priors)
samples = (v = [rand(prior) for _ in 1:10^5], weight = ones(10^5))
```

## Quantities derived from the samples

Any quantity that follows from the expected counts can be computed for each sample. Its
values over the samples are its distribution: the posterior with the samples of a fit, or
the prediction from the priors. [`get_predicted_counts`](@ref) gives the expected counts
in energy ranges for each sample, in total and for each source, and a quantity is then a
function of them. The ends of the ranges must be edges of the bins of the model. For
example, the event rate between 500 and 3000 keV, and the counts in the ``^{208}``Tl peak
at 2615 keV:

```julia
using IntervalSets, Statistics
counts = get_predicted_counts(samples, model, bin_edges, 500 .. 3000)
rate = Dict(name => uconvert.(u"s^-1", c ./ livetime) for (name, c) in counts)
peak = get_predicted_counts(samples, model, bin_edges, 2610 .. 2620)
quantile(peak[:total], [0.16, 0.5, 0.84])  # median and central 68% interval
```

[`background_index`](@ref) gives the background index at the Q-value of the 0νββ decay:
the counts in a window around it, divided by the width of the window and by the exposure,
the product of the live time and of `mass`, the mass of the HPGe detectors of the
simulation. The result is in counts/(keV kg yr). The bins of the model must cover the
window, which leaves out the γ lines at 2104 and 2119 keV. A model with one bin for each
part of the window is enough:

```julia
window_bins = [[1930, 2099], [2109, 2114], [2124, 2190]]
model_qbb, bin_edges_qbb = build_forward_model(components; norm, binning = window_bins)
bi = background_index(samples, model_qbb, bin_edges_qbb; livetime, mass = 1000u"kg")
quantile(bi[:total], 0.9)  # upper limit at 90% probability
```

[`group_histograms`](@ref) sums the values of the sources into groups, sample by sample,
for example by decay chain. The total of each group, for example the whole ``^{238}``U
background, is then a quantity like the others:

```julia
by_chain = group_histograms(bi, ["²³⁸U" => r"Pb214_to_Po214$", "²³²Th" => r"Pb212_to_Pb208$"])
quantile(by_chain[Symbol("²³⁸U")], [0.16, 0.5, 0.84])
```

## Plots

[`plot_posteriors!`](@ref) draws the distribution of a quantity for several sources or
groups, one violin per row, labelled by the names of the sources or groups or by
`labels`. The violins show the logarithm of the values, which often span orders of
magnitude. [`plot_posterior!`](@ref) draws the distribution of one of them as a
histogram, with its smallest 68% interval:

```julia
import CairoMakie
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1, 1]; xlabel = "Background index [counts/(keV kg yr)]")
plot_posteriors!(ax, by_chain)  # rows: total, ²³⁸U and ²³²Th
ax_u = CairoMakie.Axis(fig[2, 1]; xlabel = "²³⁸U background index [counts/(keV kg yr)]")
plot_posterior!(ax_u, by_chain[Symbol("²³⁸U")])
```

[`plot_spectrum!`](@ref) draws the median of the expected counts in each bin, with the
band of the central 68% interval, for the total in black and for each source or group.
With `residuals = false`, [`fit_axes`](@ref) makes a figure without the panel of the
residuals, which needs data:

```julia
fig, ax, _ = fit_axes(; residuals = false)
plot_spectrum!(ax, samples, model, bin_edges; groups = ["²³⁸U" => r"Pb214", "²³²Th" => r"Pb212"])
CairoMakie.axislegend(ax)
```
