# Accounting for systematic uncertainties

`Refit.jl` contains functionality for handling the effect of systematic uncertainties in
simulated spectra.

In a standard template fit of simulated spectra, the activity of each source scales a spectrum of fixed
shape but no systematic uncertainities on the MC shape are considered.
`Refit.jl` can incorportate systematic uncertainities in the Monte-Carlo templates directly
into the model allowing them to become fit parameters and the uncertainity accounted for.

## Generalised histograms
This is done with a [`GeneralisedHistogram`](@ref) object. This is like a `Histogram` but one
where the weights depend on parameters values.

Inside, it holds the histogram simulated at every point of a grid in the parameters. In
each bin it interpolates linearly in each parameter between the grid points around the
given values.

```@setup systematics
using Refit
using StatsBase
using LegendMakie
import Format
import CairoMakie

energy = 0.5:1.0:99.5  # bin centres in keV

# a peak on a falling continuum, x moves counts from the peak into the continuum
function spectrum(x)
    peak = exp.(-(energy .- 60) .^ 2 ./ 8)
    continuum = exp.(-energy ./ 30)
    counts = 1000 .* ((1 - x) .* peak ./ sum(peak) .+ x .* continuum ./ sum(continuum))
    return Histogram((0.0:1.0:100.0,), counts)
end

x_grid = 0.0:0.5:1.0
ghist = GeneralisedHistogram([spectrum(x) for x in x_grid], :x => x_grid)

fig = CairoMakie.Figure(size = (1100, 450))
ax = CairoMakie.Axis(fig[1, 1]; xlabel = "Energy [keV]", ylabel = "Counts per keV")
for x in (0.0, 0.5, 1.0)
    plot_hist!(ax, get_histogram(ghist, (x = x,)); label = "x = $x, grid point")
end
plot_hist!(
    ax,
    get_histogram(ghist, (x = 0.25,));
    color = :black,
    linestyle = :dash,
    label = "x = 0.25, interpolated",
)
CairoMakie.Legend(fig[1, 2], ax; framevisible = false)
```

As an example, take a spectrum with a peak on a falling continuum, where a parameter `x`
moves counts from the peak into the continuum. Its histogram is simulated at three grid
points, `x` = 0, 0.5 and 1. The generalised histogram `ghist` gives it at any `x` in
between: `get_histogram(ghist, (x = 0.25,))` is the dashed line.

```@example systematics
fig # hide
```

For given values of the parameters, [`get_histogram`](@ref) returns a `Histogram` and
[`get_weights`](@ref) its counts per bin. The values are given as a named tuple, in any
order. Other entries, such as the activities of a fit, are ignored. A value outside the
grid is an error.

## Reading a generalised histogram

A generalised histogram needs:

- a grid of evenly spaced values for each parameter, with any number of parameters,
- the histogram simulated at every point of the grid, all with the same binning.

The histograms usually come from a simulation that writes them to an LH5 file, in the
format described [below](@ref "File format of a generalised histogram").
[`read_generalised_hist`](@ref) reads one. [`read_generalised_pdf`](@ref) also divides
by the number of simulated primaries and averages several files with weights, like
[`read_pdf`](@ref):

```julia
ghist = read_generalised_pdf("simulation.lh5"; groups = "pdf/mul1/all")
```

!!! note "Building a generalised histogram from histograms"
    A generalised histogram can also be built directly from histograms in Julia, as in
    the example above. `hists` is an array with one dimension per parameter, and each
    further argument is a pair `name => values` for one dimension, in the order of the
    dimensions:

    ```julia
    ghist = GeneralisedHistogram(hists, :x => 0.0:0.5:1.0)
    ```

## File format of a generalised histogram

A generalised histogram is stored in an LH5 group with two parts: the description of
the grid in `grid_info`, and the histogram of every grid point in `grid`. In a file of
simulated spectra it takes the place of a histogram, for example at `pdf/mul1/all`,
next to the number of simulated primaries. For two parameters, whose names stand here
as `<name0>` and `<name1>`:

```text
/
├── nr_sim_events      number of simulated primaries
└── pdf
    └── mul1
        └── all
            ├── grid_info
            │   ├── axis_0
            │   │   ├── name      "<name0>"
            │   │   └── values
            │   │       ├── first
            │   │       ├── last
            │   │       └── step
            │   └── axis_1
            │       ├── name      "<name1>"
            │       └── values    first, last, step
            └── grid
                ├── <name0>_0
                │   ├── <name1>_0     histogram
                │   ├── <name1>_1     histogram
                │   :
                ├── <name0>_1
                :
```

The names are free, for example `fccd` and `tl` in the [FCCD tutorial](tutorial_fccd.md).

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

[`get_max_activity`](@ref) also takes a generalised histogram, to set the range of its
activity from another energy window by hand. It returns the largest activity over the
grid points. Here a made-up measured spectrum, the simulation at `x = 0.3` with an
activity of 2, sets the range from the counts of the peak between 50 and 70 keV:

```@example systematics
using IntervalSets

data = Histogram((0.0:1.0:100.0,), 2 .* get_weights(ghist, (x = 0.3,)))
priors = Dict(:source => 0.0 .. 1.5 * get_max_activity(50, 70, data, ghist))
```

The grid point `x = 1`, with the fewest counts in the peak, gives the largest activity.

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
