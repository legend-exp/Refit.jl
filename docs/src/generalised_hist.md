# Generalised histograms

## Motivation

In a fit of simulated spectra, the activity of each source scales a spectrum of fixed
shape. Some properties of the setup change the shape itself, for example the thickness
of the dead layer of a detector (FCCD, full charge collection depth), the width of its
transition layer, or the position of a source. A scale factor cannot describe them.

The way around this is to simulate the spectrum on a grid of values of these
properties and to interpolate between the grid points, also called template morphing.
The properties then become parameters of the fit, like the activities.

## What it does

A [`GeneralisedHistogram`](@ref) holds the spectrum simulated at every point of an evenly
spaced grid in 1 to 5 parameters. [`get_weights`](@ref) returns the counts per bin at any
point inside the grid. In every bin it interpolates linearly in each parameter between
the neighbouring grid points:

```julia
gh = read_generalised_pdf("templates.lh5"; groups = "pdf/geds")
get_weights(gh, (fccd = 1.05, tl = 0.27))  # counts per bin, here per simulated primary
```

The user only gives the values of the parameters. [`rebin_hist`](@ref) and
[`add_histograms`](@ref) work as for a `Histogram`, on every grid point at once. A
generalised histogram is read from an LH5 file with [`read_generalised_hist`](@ref), or
with [`read_generalised_pdf`](@ref), which also divides by the number of simulated
primaries like [`read_pdf`](@ref). The docstring of `read_generalised_hist` describes the
format of the file.

The [FCCD tutorial](tutorial_fccd.md) fits the FCCD and the transition layer of a
detector in this way.

## Limitations

- The interpolation is linear in each bin. This suits a parameter that changes the
  counts in each bin smoothly, such as a dead layer or a source position. It does not
  suit a parameter that moves a peak in energy, such as the energy scale: between two
  grid points the peak would appear twice, at half the height. Such parameters belong in
  an [analytic part](@ref "Analytic functions") of the model.
- Between grid points the interpolation is an approximation, unless the counts in each
  bin are linear in the parameter. Check that the grid is fine enough, for example by
  comparing with a simulation halfway between two grid points.
- A point outside the grid is an error. The prior range of each parameter must lie
  inside the grid.
