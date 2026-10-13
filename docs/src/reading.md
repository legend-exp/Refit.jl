# Reading simulated spectra

The simulated spectra are stored in LH5 files. `read_hist(group, file)` reads one
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
name of the field with the number of primaries.

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
