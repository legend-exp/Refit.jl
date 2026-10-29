# Expected background from known activities

Without a measured spectrum, the activities of the components of an experiment can come
from radioassay measurements, upper limits or decay rates. Drawing activities from their
distributions and passing each draw through the [forward model](forward_model.md) gives
the distribution of the expected spectrum and of quantities derived from it, such as the
background index at the Q-value of the 0νββ decay. The same functions take the samples
of a fit, once there is a measured spectrum.

## Activities

[`activity_prior`](@ref) turns a measured activity, with its uncertainty, into a Gaussian
cut at 0, and an upper limit into a half-normal distribution whose quantile at the
confidence level is the limit:

```julia
using Unitful
activity_prior(; val = 83u"μBq", unc = 113u"μBq")
activity_prior(; limit = 12u"μBq", cl = 0.9)
```

[`read_activities`](@ref) reads them from a YAML file, which lists the components of the
experiment, their quantity (mass, area or number of pieces) and the activity per unit of
quantity of each simulation (`simid`) that stands for them. The activities are the
quantity times the activity per unit, in µBq. A production cycle of legend-simflow has
the file in its configuration directory:

```julia
prod = SimProd("prod/l1000dsg01/v1.1.0")
activities = read_activities(prod)  # Dict(component => Dict(simid => distribution))
```

The result is a dictionary, so components can be added or changed by hand, for example
with an activity from a half-life:

```julia
# 2νββ decay of the 7.3e27 nuclei of 76Ge in 1000 kg of germanium enriched to 92%, with a
# half-life of 2.02e21 yr, with an uncertainty of 2%
a = log(2) / 2.02e21u"yr" * 7.3e27
activities[:two_nu_bb] = Dict(:hpge_icpc_bulk_2vbb => activity_prior(; val = a, unc = 0.02 * a))
```

## Sources and priors

[`activity_components`](@ref) makes one source of the forward model for each component
and simid, named `<component>_<simid>`, with the spectrum of the simid and the
distribution of the activity of the component. A simulation that stands for parts of
different materials gives one source for each part. With the production cycle, the
spectra are read with the event selection of [`read_pdf`](@ref):

```julia
components, priors = activity_components(prod, activities; groups = [:is_ss_aoe, :not_spms_coincident])
```

`norm` turns an activity in µBq into the number of decays in the live time. With the
live time as a quantity with units (from Unitful.jl), it is the number µBq × live time:

```julia
livetime = 10u"yr"
norm = ustrip(NoUnits, u"μBq" * livetime)
model, bin_edges = build_forward_model(components; norm, e_low = 500, e_high = 3000, binning = 10)
```

The activities are drawn from the priors with BAT.jl. There is no likelihood, so the
draws are independent:

```julia
using BAT
samples = bat_sample(build_prior_dist(priors), IIDSampling(nsamples = 10^5)).result
```

## Background index and other quantities

[`get_predicted_counts`](@ref) gives the counts predicted by the model in energy ranges,
for each draw, in total and for each source. The ends of the ranges must be edges of the
bins of the model. Any quantity follows from the counts, for example the event rate:

```julia
using IntervalSets
counts = get_predicted_counts(samples, model, bin_edges, 500 .. 3000)
rate = Dict(name => uconvert.(u"s^-1", c ./ livetime) for (name, c) in counts)
```

[`background_index`](@ref) divides the counts in the window around the Q-value of the
0νββ decay by the width of the window and by the exposure, the product of the live time
and of `mass`, the mass of the HPGe detectors of the simulation. The result is in
counts/(keV kg yr). The bins of the model must cover the window, which leaves out the γ
lines at 2104 and 2119 keV. A model with one bin for each part of the window is enough:

```julia
using Statistics
window_bins = [[1930, 2099], [2109, 2114], [2124, 2190]]
model_qbb, bin_edges_qbb = build_forward_model(components; norm, binning = window_bins)
bi = background_index(samples, model_qbb, bin_edges_qbb; livetime, mass = 1000u"kg")
quantile(bi[:total], [0.16, 0.5, 0.84])  # median and central 68% interval
quantile(bi[:total], 0.9)                # upper limit at 90% probability
```

[`group_histograms`](@ref) sums the sources into groups, draw by draw, and
[`plot_posteriors!`](@ref) draws the distribution of each source or group as a violin, one
row each. The violins show the logarithm of the values, which span orders of magnitude:

```julia
import CairoMakie
by_chain = group_histograms(bi, ["²¹⁴Pb" => r"Pb214_to_Po214$", "²¹²Pb" => r"Pb212_to_Pb208$"])
fig = CairoMakie.Figure()
ax = CairoMakie.Axis(fig[1, 1]; xlabel = "Background index [counts/(keV kg yr)]")
plot_posteriors!(ax, by_chain)
```

## Expected spectrum

[`plot_spectrum!`](@ref) draws the median of the expected counts in each bin, with the
band of the central 68% interval, for the total in black and for each source or group.
With `residuals = false`, [`fit_axes`](@ref) makes a figure without the panel of the
residuals, which needs data:

```julia
fig, ax, _ = fit_axes(; residuals = false)
plot_spectrum!(ax, samples, model, bin_edges; groups = ["²¹⁴Pb" => r"Pb214", "²¹²Pb" => r"Pb212"])
CairoMakie.axislegend(ax)
```
