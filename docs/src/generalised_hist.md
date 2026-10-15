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
Inside, it holds the histogram simulated at every point of a grid in the parameters. In
each bin it interpolates linearly in each parameter between the grid points around the
given values.

For example, the [FCCD tutorial](tutorial_fccd.md) uses the spectra of a 60 keV γ line,
simulated on a grid of FCCD and transition layer thickness (TL), in mm.
[`read_generalised_pdf`](@ref) reads them from a file and divides them by the number of
simulated γs:

```@example systematics
using Refit
using StatsBase

file = "../../test/test_files/fccd/templates.lh5"  # relative to docs/src
ghist = read_generalised_pdf(file; groups = "pdf/geds")
ghist.grid
```

For given values of the parameters, [`get_histogram`](@ref) returns a `Histogram` and
[`get_weights`](@ref) its counts per bin. The values are given as a named tuple, in any
order. Other entries, such as the activities of a fit, are ignored. A value outside the
grid is an error. The spectra for three thicknesses of the transition layer, which sets
the level of the continuum below the peak, are below. The thickness of 0.25 mm lies
between two grid points:

```@example systematics
using LegendMakie
import Format
import CairoMakie

fig = CairoMakie.Figure(size = (1100, 450))
ax = CairoMakie.Axis(
    fig[1, 1];
    yscale = log10,
    xlabel = "Energy [keV]",
    ylabel = "Counts per γ and keV",
)
for tl in (0.1, 0.25, 0.4)  # mm
    h = get_histogram(ghist, (fccd = 1.0, tl = tl))
    plot_hist!(ax, h; label = "FCCD = 1.0 mm, TL = $tl mm")
end
CairoMakie.xlims!(ax, 0, 100)
CairoMakie.ylims!(ax, 1e-4, 1)
CairoMakie.Legend(fig[1, 2], ax; framevisible = false)
fig
```

## Building a generalised histogram

A generalised histogram needs:

- a grid of evenly spaced values for each parameter, with any number of parameters,
- the histogram simulated at every point of the grid, all with the same binning.

The histograms usually come from a simulation that writes them to an LH5 file, in the
format described [below](@ref "File format of a generalised histogram").
[`read_generalised_hist`](@ref) reads one. [`read_generalised_pdf`](@ref) also divides
by the number of simulated primaries and averages several files with weights, like
[`read_pdf`](@ref).

For custom cases, a generalised histogram can also be built in Julia from an array of
histograms, with one dimension per parameter. Each further argument is a pair
`name => values` for one dimension, in the order of the dimensions. For example, the
spectrum of a 60 keV γ line for a dead layer of thickness `d` in mm, without a
transition layer, if the γs interact uniformly in the first 2 mm of the detector: a
fraction (2 − d) / 2 is in the peak, the rest is lost.

```@example systematics
function line_spectrum(d)
    weights = zeros(100)  # 1 keV bins from 0 to 100 keV
    weights[61] = (2 - d) / 2  # the bin from 60 to 61 keV
    return Histogram((0.0:1.0:100.0,), weights)
end

dead_layer = 0.5:0.25:1.5  # mm
custom = GeneralisedHistogram([line_spectrum(d) for d in dead_layer], :d => dead_layer)
get_weights(custom, (d = 0.6,))[61]
```

## File format of a generalised histogram

A generalised histogram is stored in an LH5 group with two parts: the description of
the grid in `grid_info`, and the histogram of every grid point in `grid`. In a file of
simulated spectra it takes the place of a histogram, for example at `pdf/geds`, next to
the number of simulated primaries:

```text
/
├── nr_sim_events      number of simulated primaries
└── pdf
    └── geds
        ├── grid_info
        │   ├── axis_0
        │   │   ├── name      "fccd"
        │   │   └── values
        │   │       ├── first     0.6
        │   │       ├── last      1.4
        │   │       └── step      0.2
        │   └── axis_1
        │       ├── name      "tl"
        │       └── values    first, last, step
        └── grid
            ├── fccd_0
            │   ├── tl_0      histogram
            │   ├── tl_1      histogram
            │   :
            ├── fccd_1
            :
```

The conventions are:

- Each `axis_d` of `grid_info`, counted from 0, is one parameter: its `name` and evenly
  spaced `values` from `first` to `last` in steps of `step`. The values may carry their
  units as an attribute, which is not read.
- The groups in `grid` are nested in the order of the axes. The group `<name>_<i>` holds
  the `i`-th value of that parameter, counted from 0. Every combination of values must
  be present.
- Each grid point is an LH5 histogram. All of them have the same binning and hold the
  raw counts of the simulation.
- Every grid point is simulated with the same number of primaries, `nr_sim_events`.
  This is required: a single number normalises all the grid points.

## Selecting the fit range

[`select_bins`](@ref) selects the bins of a fit in the same way as
[`build_forward_model`](@ref): it rebins with `binning` (see [Binning](binning.md)) and
keeps the bins fully inside `e_low` to `e_high`, at every grid point. It works the same
for a `Histogram`, such as the measured spectrum, so both get the same bins:

```@example systematics
ghist_fit, bin_edges = select_bins(ghist; binning = 5, e_low = 5, e_high = 100)
get_weights(ghist_fit, (fccd = 1.05, tl = 0.27))
```

Selecting bins and interpolating both sum counts, so the order does not matter. Select
the bins once, before the fit, so that each evaluation only touches the fitted bins.
[`rebin_hist`](@ref) rebins every grid point without selecting a range.

## Fits with a generalised histogram

[`build_forward_model`](@ref) takes a generalised histogram in place of a `Histogram` in
any component (see [Model components](@ref)). At each evaluation it interpolates the
spectrum at the values of the parameters, then scales it by the activity:

```julia
components = Dict(:source => ghist, :background => (analytic = (E, p) -> p.bkg,))
model, bin_edges = build_forward_model(components; e_low = 5, e_high = 100, binning = 5)
```

[`get_priors`](@ref) estimates the range of the activity at every grid point and uses
the widest. The prior range of each parameter of the grid must lie inside its grid.

For a simple fit, the model can also be written by hand with the bins of
[`select_bins`](@ref), as in the [FCCD tutorial](tutorial_fccd.md). Any function of the
parameters that returns the expected counts in the bins of the measured counts is a
valid forward model for [`build_likelihood`](@ref).

!!! warning "The interpolation is an approximation"
    Between grid points, linear interpolation is exact only if the counts in each bin are
    linear in the parameter. It suits a parameter that changes the counts in each bin
    smoothly, such as a dead layer or a source position. It does not suit a parameter
    that moves a peak in energy, such as the energy scale: between two grid points the
    peak appears twice, at half the height. Such parameters belong in an
    [analytic part](@ref "Analytic functions") of the model. It is up to the user to
    check that the grid is fine enough, for example by comparing with a simulation
    halfway between two grid points.
