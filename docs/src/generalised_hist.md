# Generalised histograms

Some simulated spectra depend on a parameter that is not known well enough to fix, for
example the thickness of a dead layer or the position of a source. A
[`GeneralisedHistogram`](@ref) holds the spectrum simulated at every point of a grid in
up to 5 such parameters. It returns the spectrum at any point inside the grid by linear
interpolation between the grid points, separately in each bin.

## Building a generalised histogram

The spectra are given as an array of histograms with one dimension per parameter, all
with the same binning. The keyword arguments give the name and the grid values of each
parameter, in the order of the dimensions. Each grid is an evenly spaced range.

As a toy example, a flat spectrum from 0 to 100 keV with a peak whose height falls with
the FCCD (in mm), simulated at 4 FCCD values times 3 values of a second parameter
`alpha_frac` that scales the flat part:

```@example ghist
using Refit, StatsBase

fccd = 0.5:0.5:2.0
alpha_frac = 0.0:0.5:1.0

function toy_spectrum(f, a)
    weights = fill(10.0 * (1 + a), 100)
    weights[50] += 1000 / f
    return Histogram((0.0:1.0:100.0,), weights)
end

hists = [toy_spectrum(f, a) for f in fccd, a in alpha_frac]
gh = GeneralisedHistogram(hists; fccd = fccd, alpha_frac = alpha_frac)
nothing # hide
```

`hists[i, j]` is the spectrum at the `i`-th FCCD value and the `j`-th value of
`alpha_frac`.

## Spectrum at a point of the grid

[`get_weights`](@ref) returns the counts per bin at any point inside the grid, and
[`get_histogram`](@ref) the same as a `Histogram`. The parameters are given as a named
tuple, in any order. Other entries, such as the activities of a fit, are ignored:

```@example ghist
w = get_weights(gh, (fccd = 0.75, alpha_frac = 0.25, A = 3.0))
w[49:51]
```

Halfway between the FCCD values 0.5 and 1.0 mm, the peak is halfway between 2000 and
1000 counts. A point outside the grid is an error.

## Rebinning

[`rebin_hist`](@ref) rebins the spectrum at every grid point, in the same ways as for a
`Histogram` (see [Binning](binning.md)). Rebinning to edges keeps only the range of the
edges, so it also selects the fitted energy range:

```@example ghist
gh_fit = rebin_hist(gh, collect(40.0:5.0:60.0))
get_weights(gh_fit, (fccd = 0.75, alpha_frac = 0.25))
```

Rebinning and interpolation both sum counts, so the order does not matter. Rebin once,
before the fit, so that each evaluation touches only the fitted bins.

## Limitations

The interpolation is linear in each bin. This suits parameters that change the counts in
each bin smoothly, such as a dead layer or a source position. It does not suit a
parameter that moves a peak in energy, such as the energy scale: between two grid points
the peak would appear twice, at half the height. Such parameters belong in an analytic
part of the model instead.
