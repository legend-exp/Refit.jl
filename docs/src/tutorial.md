# Tutorial: fitting a background spectrum

This tutorial fits a measured energy spectrum with a sum of simulated spectra, from
reading the files to sampling the posterior. It runs every time the documentation
is built.

The files are made up. A script, `docs/src/tutorial_data.jl`, writes them in LH5
format to a temporary directory. They stand in for the output of a real simulation
and for a measured spectrum. Each source with a gamma line has a peak with a 1 keV
width on a flat or linear continuum below the peak. The 2νββ spectrum is a broad
Gaussian (600 keV width, centered at 1 MeV).

- one file for each simulated source with the histograms `pdf/mul1/icpc` and
  `pdf/mul1/bege` (the two detector types) and the number of simulated primaries,
- one file `data.lh5` with the measured histograms `hist/mul1/icpc` and
  `hist/mul1/bege`.

The measured spectrum is generated from known activities, so the result of the fit
can be compared with the truth at the end.

```@setup tutorial
include("tutorial_data.jl")
```

## Settings

```@example tutorial
using Refit
using StatsBase

livetime = 1.0e6  # s
e_low, e_high = 100, 4000  # keV, the range of the fit
binning = 10  # merge 10 bins of 1 keV into one
blind = (2015.0, 2065.0)  # keV, a region around the 2039 keV Q-value that is not fitted

# The activities are in mBq and the simulated spectra are counts per simulated
# primary, so the expected counts are 1e-3 * livetime * activity * spectrum
norm = 1e-3 * livetime

# the histograms of the measured spectrum and of each simulation are summed over these groups
spec = ["mul1/icpc", "mul1/bege"]
nothing # hide
```

## Reading the data

`read_hist` reads one histogram. The data of the two detector types are added with
`add_histograms`.

```@example tutorial
data = add_histograms([read_hist("hist/$s", joinpath(dir, "data.lh5")) for s in spec])
```

## Reading the simulated spectra

Each source has its own file. `read_pdf` sums the groups and divides by the number of
simulated primaries, so that the spectrum is in counts per primary. The name of each
spectrum is a `Symbol`, the activity of the source gets the same name.

```@example tutorial
groups = ["pdf/$s" for s in spec]

components = ["cables", "fibers"]
isotopes = ["K40", "Pb212_to_Pb208"]

mc_dict = Dict(
    Symbol(comp, "_", iso) => read_pdf(joinpath(dir, "$(comp)_$iso.lh5"); groups) for
    comp in components for iso in isotopes
)
keys(mc_dict)
```

The 2νββ decay is simulated once for each detector type. `read_pdf` accepts a
dictionary from file name to weight and returns the weighted average of the spectra.
The weights are normalised. Here the weight is the mass of each detector type in kg.

```@example tutorial
mass = Dict("icpc" => 80.0, "bege" => 20.0)
bulk = Dict(joinpath(dir, "hpge_$(t)_bulk_2vbb.lh5") => w for (t, w) in mass)

mc_dict[:two_vbb_hpge_bulk] = read_pdf(bulk; groups)
nothing # hide
```

Alpha events have no simulation. As a simple model they are a constant spectrum. Its
activity is the number of counts in each 1 keV bin, so the spectrum is
`1 / norm` in every bin.

```@example tutorial
mc_dict[:alpha] = Histogram(data.edges, fill(1 / norm, length(data.weights)))
nothing # hide
```

## Forward model and likelihood

The [forward model](forward_model.md) maps the activities to the expected counts in
each energy bin. The [likelihood](likelihood.md) compares these with the measured
counts, which are taken from the data in the same bins.

```@example tutorial
model, bin_edges = build_forward_model(mc_dict; e_low, e_high, blind, binning, norm)
counts = get_counts(data, bin_edges)
likelihood = build_likelihood(counts, model)
nothing # hide
```

## Ranges for the activities

The fit needs a range for each activity. `get_priors` estimates it from the data
(see [Priors](priors.md)). The reference ranges give the energy window of a
peak of each isotope in keV. The 1461 keV line of ``^{40}``K and the 2615 keV line of
``^{208}``Tl limit the activities of these sources.

```@example tutorial
ref_ranges = Dict(:K40 => (1456, 1466), :Pb212_to_Pb208 => (2610, 2620))

ranges = get_priors(mc_dict, data; ref_ranges, e_low, e_high, blind, norm)
```

A range can be changed by hand where the guess is poor, for example
`ranges[:two_vbb_hpge_bulk] = 0.0 .. 20.0` (with `using IntervalSets`).

The prior of each activity is a uniform distribution on its range. The prior is a
`NamedTupleDist` from ValueShapes.jl, with the same names as the activities.

```@example tutorial
using Distributions
using ValueShapes

prior = NamedTupleDist(; (k => Uniform(r.left, r.right) for (k, r) in ranges)...)
nothing # hide
```

## Sampling the posterior

[BAT.jl](https://bat.github.io/) samples the posterior with a Markov chain Monte
Carlo (a random walk through the space of the parameters that visits each value
with a probability given by the posterior). The first steps of each chain are
discarded as burn-in. This run is short, a real fit uses more steps.

```@example tutorial
using BAT

posterior = PosteriorMeasure(likelihood, prior)

samples = bat_sample(
    posterior,
    TransformedMCMC(
        proposal = RandomWalk(),
        strict = false,
        burnin = MCMCMultiCycleBurnin(nsteps_per_cycle = 10_000, max_ncycles = 30),
        nsteps = 10^5,
        nchains = 4,
    ),
).result
nothing # hide
```

## Result

The most probable value of each activity is compared with the one used to make the
measured spectrum.

```@example tutorial
best = mode(samples)
for name in keys(true_activities)
    println(rpad(name, 25), "fit ", round(best[name], digits = 1),
            " mBq, true ", true_activities[name], " mBq")
end
println(rpad(:alpha, 25), "fit ", round(best.alpha, digits = 2),
        " counts per keV, true ", true_alpha)
```
