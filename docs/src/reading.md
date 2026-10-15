# Reading simulated spectra

The simulated and data spectra are read from LH5 files (see the [LH5 spec](https://legend-exp.github.io/legend-data-format-specs/dev/hdf5/)). `read_hist(group, file)` reads one
histogram from the group `group` of `file`.

The examples below use a small file with two histograms and the number of simulated
primaries (`nr_sim_events`), written here from a test file:

```@setup demo
using Refit, LegendHDF5IO

source = "../../test/test_files/test_histograms.lh5"
binning, weights = lh5open(source, "r") do f
    f["hist/binning/axis_0"], f["hist/weights"][:]
end

file = tempname() * ".lh5"
lh5open(file, "w") do f
    f["pdf/mul_surv/binning/axis_0"] = binning
    f["pdf/mul_surv/weights"] = weights
    f["pdf/mul_other/binning/axis_0"] = binning
    f["pdf/mul_other/weights"] = 2 .* weights
    f["nr_sim_events"] = 1000
end
```

```@example demo
# one raw histogram from the group "pdf/mul_surv"
h = read_hist("pdf/mul_surv", file)
```

## Counts per simulated primary

A simulated spectrum is normally used as counts per simulated primary, so that the
activity of the source scales it to the measured counts (see the
[forward model](forward_model.md)). `read_pdf` reads the histogram and divides it by
the number of primaries in the field `nr_sim_events`:

```@example demo
# the default group, divided by the number of primaries
pdf = read_pdf(file)
```

`groups` is the name of the group to read, `"pdf/mul_surv"` by default. A vector of
names reads several groups, which are summed before the division. `n_prim_field` is the
name of the field with the number of primaries. With `n_prim_field = nothing` the
histogram is not divided, and the result is the same as from `read_hist`.

```@example demo
# sum two groups, then divide by the number of primaries
pdf = read_pdf(file; groups = ["pdf/mul_surv", "pdf/mul_other"])
```

## Weighting the groups

`weights` gives one factor for each group, for example the fraction of the decays in
each component. The groups are multiplied by their weight before they are summed:

```@example demo
pdf = read_pdf(file; groups = ["pdf/mul_surv", "pdf/mul_other"], weights = [0.7, 0.3])
```

Any other factor can be applied to the result, which is a normal `Histogram`. For
example, to divide by 2 in place:

```julia
pdf.weights ./= 2
```

## The same source in several files

A source is often simulated once per detector type. `read_pdf` takes a dictionary from
file name to weight, for example the mass of each detector type, and returns the
weighted average of the spectra. The weights are divided by their sum, and each file
is divided by its own number of primaries.

```julia
files = Dict(
    "icpc_bulk_2vbb-tier_pdf.lh5" => 87.1,  # detector mass in kg
    "bege_bulk_2vbb-tier_pdf.lh5" => 19.1,
)
pdf = read_pdf(files; groups = "pdf/mul_surv")
```
