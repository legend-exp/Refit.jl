# Accounting for systematic uncertainties

`Refit.jl` contains functionality for handling the effect of systematic uncertainties in
simulated spectra.

In a fit of simulated spectra, the activity of each source scales a spectrum of fixed
shape. Some properties of the setup are only known within an uncertainty and change the
shape itself, for example the thickness of the dead layer of a detector (FCCD, full
charge collection depth), the width of its transition layer, or the position of a
source. A scale factor cannot describe them. Instead, the spectrum is simulated on a
grid of values of these properties and interpolated between the grid points, also
called template morphing. The properties then become parameters of the fit, and their
uncertainty is part of the result.

## Generalised histograms

A [`GeneralisedHistogram`](@ref) is like a `Histogram` whose counts depend on parameters.
For given values of the parameters, [`get_histogram`](@ref) returns a `Histogram` and
[`get_weights`](@ref) its counts per bin:

```julia
get_histogram(ghist, (fccd = 1.05, tl = 0.27))  # a Histogram
get_weights(ghist, (fccd = 1.05, tl = 0.27))    # its counts per bin
```

Inside, it holds the histogram simulated at every point of a grid in the parameters. In
each bin it interpolates linearly in each parameter between the grid points around the
given values.

## Building a generalised histogram

A generalised histogram needs:

- a grid of evenly spaced values for each parameter, with any number of parameters,
- the histogram simulated at every point of the grid, all with the same binning.

The histograms usually come from a simulation that writes them to an LH5 file.
[`read_generalised_hist`](@ref) reads one, and its docstring describes the format of the
file. [`read_generalised_pdf`](@ref) also divides by the number of simulated primaries,
like [`read_pdf`](@ref). For example, the simulated spectra of the
[FCCD tutorial](tutorial_fccd.md) are the spectra of a 60 keV γ line on a grid of FCCD
and transition layer thickness (TL), in mm:

```@example systematics
using Refit
using StatsBase

file = "../../test/test_files/fccd/templates.lh5"  # relative to docs/src
ghist = read_generalised_pdf(file; groups = "pdf/geds")
ghist.grid
```

A generalised histogram can also be built in Julia, from a function that returns the
histogram for given values of the parameters. The keyword arguments give the grid of
each parameter. As a toy example, a flat spectrum with a peak whose height falls with
the FCCD:

```@example systematics
function toy_spectrum(p)
    weights = fill(10.0, 100)  # 1 keV bins from 0 to 100 keV
    weights[50] += 1000 / p.fccd
    return Histogram((0.0:1.0:100.0,), weights)
end

toy = GeneralisedHistogram(toy_spectrum; fccd = 0.5:0.5:2.0)
nothing # hide
```

## Spectrum at given parameter values

The values of the parameters are given as a named tuple, in any order. Other entries,
such as the activities of a fit, are ignored. A value outside the grid is an error.

```@example systematics
w = get_weights(ghist, (tl = 0.27, fccd = 1.05, activity = 3.0))
w[58:62]
```

## Rebinning and the fit range

[`rebin_hist`](@ref) rebins the histogram at every grid point, in the same ways as a
`Histogram` (see [Binning](binning.md)). Rebinning to bin edges keeps only the range of
the edges, so it also selects the fit range:

```@example systematics
ghist_fit = rebin_hist(ghist, collect(5.0:5.0:100.0))  # 5 to 100 keV in 5 keV bins
get_weights(ghist_fit, (fccd = 1.05, tl = 0.27))
```

Rebinning and interpolation both sum counts, so the order does not matter. Rebin once,
before the fit, so that each evaluation only touches the fitted bins.

## Fits with a generalised histogram

The model of a fit is a function of the parameters that returns the expected counts in
each bin. With a generalised histogram it interpolates at the parameter values of each
evaluation, for example for a known number of γs:

```julia
model(p) = n_gamma .* get_weights(ghist_fit, p)
likelihood = build_likelihood(counts, model)
```

[`build_likelihood`](@ref) takes any such function (see [Likelihood](likelihood.md)).
The prior range of each parameter must lie inside its grid. The
[FCCD tutorial](tutorial_fccd.md) is a complete example.

!!! warning "The interpolation is an approximation"
    Between grid points, linear interpolation is exact only if the counts in each bin are
    linear in the parameter. It suits a parameter that changes the counts in each bin
    smoothly, such as a dead layer or a source position. It does not suit a parameter
    that moves a peak in energy, such as the energy scale: between two grid points the
    peak appears twice, at half the height. Such parameters belong in an
    [analytic part](@ref "Analytic functions") of the model. It is up to the user to
    check that the grid is fine enough, for example by comparing with a simulation
    halfway between two grid points.
