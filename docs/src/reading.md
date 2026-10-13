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
h = read_hist("pdf/mul_surv", file)
```

## Counts per simulated primary

A simulated spectrum is normally used as counts per simulated primary, so that the
activity of the source scales it to the measured counts (see the
[forward model](forward_model.md)). `read_pdf` reads the histogram and divides it by
the number of primaries in the field `nr_sim_events`:

```@example demo
pdf = read_pdf(file)
```

The histogram is counts per primary, here `h.weights ./ 1000`:

```@example demo
pdf.weights ≈ h.weights ./ 1000
```

`groups` lists several groups, which are summed before the division. `n_prim_field`
is the name of the field with the number of primaries:

```@example demo
pdf = read_pdf(file; groups = ["pdf/mul_surv", "pdf/mul_other"], n_prim_field = "nr_sim_events")
pdf.weights ≈ 3 .* h.weights ./ 1000
```

## Scaling and summing by hand

The weights of a histogram can be changed directly. This is useful for any other factor,
for example a branching ratio or the fraction of decays in a component:

```@example demo
pdf = read_pdf(file)
pdf.weights ./= 2   # the histogram is divided in place
pdf.weights ≈ h.weights ./ 2000
```

`add_histograms` sums histograms with the same binning. Scaling each one before
gives a weighted sum, for example to merge two components of a detector:

```@example demo
pdf_a = read_pdf(file; groups = ["pdf/mul_surv"])
pdf_b = read_pdf(file; groups = ["pdf/mul_other"])

pdf_a.weights .*= 0.7
pdf_b.weights .*= 0.3

merged = add_histograms([pdf_a, pdf_b])
merged.weights ≈ (0.7 .* h.weights .+ 0.3 .* 2 .* h.weights) ./ 1000
```
