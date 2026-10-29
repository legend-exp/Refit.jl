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

## Spectra from a legend-simflow production cycle

[legend-simflow](https://github.com/legend-exp/legend-simflow) writes one `pdf` tier
file per simulation (`simid`) of a production cycle. [`SimProd`](@ref) finds
these files from the directory of the cycle, the one with `simflow-config.yaml`:

```julia
prod = SimProd("prod/l1000dsg01/v1.1.0")
simids(prod)  # the simulations in the cycle, as Symbols
pdf_file(prod, :hpge_icpc_bulk_2vbb)  # the path of one pdf tier file
```

All files of a cycle contain the same histograms. The `pdf` tier sorts the events into
categories: by multiplicity (`mul1`), by the LAr veto (`is_spms_coincident`,
`not_spms_coincident`) and by the A/E cut (`is_high_aoe`, `is_ss_aoe`, `is_low_aoe`,
`is_invalid_aoe`), for example. It stores one histogram for each combination of
categories and each HPGe group. [`cut_groups`](@ref) lists the combinations and
[`hpge_groups`](@ref) the HPGe groups, read from the file of the first
simulation. [`pdf_groups`](@ref) lists all the histograms of the `pdf` tier. Showing
`prod` in the REPL prints a summary of the cycle:

```julia
julia> prod
SimProd at /data/l1000dsg01/v1.1.0
├─ experiment: l1000dsg01
├─ 36 simulations:
│    cable_connectors_Pb212_to_Pb208
│    cable_connectors_Pb214_to_Po214
│    cables_Pb212_to_Pb208
│    ⋮
│    wls_reflector_Pb214_to_Po214
├─ cut groups:
│    is_high_aoe ∧ is_spms_coincident ∧ mul1
│    is_invalid_aoe ∧ is_spms_coincident ∧ mul1
│    is_low_aoe ∧ is_spms_coincident ∧ mul1
│    is_spms_coincident ∧ is_ss_aoe ∧ mul1
│    is_high_aoe ∧ mul1 ∧ not_spms_coincident
│    is_invalid_aoe ∧ mul1 ∧ not_spms_coincident
│    is_low_aoe ∧ mul1 ∧ not_spms_coincident
│    is_ss_aoe ∧ mul1 ∧ not_spms_coincident
└─ HPGe groups: all
```

The simulations, categories and HPGe groups are `Symbol`s. `read_pdf` sums the
histograms whose combination contains all the categories of `groups`, and divides by the
number of primaries. `groups` is a category or a vector of categories, in any order. A
vector of such vectors sums the histograms matched by any of them, each counted once.
`hpge_group` is the HPGe group, `:all` by default, or a vector of HPGe groups to sum.
Groups that share detectors, such as `:all` and any other, count their events twice:

```julia
# all multiplicity-1 events of the HPGe group "all"
pdf = read_pdf(prod, :hpge_icpc_bulk_2vbb)
# after the LAr veto and the A/E cuts
pdf = read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = [:not_spms_coincident, :is_ss_aoe])
# after the A/E cuts only
pdf = read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = :is_ss_aoe)
# after the low-side A/E cut only: the histograms matched by any of the lists
pdf = read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = [[:is_ss_aoe], [:is_high_aoe]])
```

[`has_pdf`](@ref) checks whether a spectrum is there without reading it: a simulation,
a selection with the same `groups` and `hpge_group` as `read_pdf`, or a histogram path:

```julia
has_pdf(prod, :hpge_icpc_bulk_2vbb)
has_pdf(prod, :hpge_icpc_bulk_2vbb; groups = [:is_ss_aoe, :not_spms_coincident], hpge_group = :all)
has_pdf(prod, :hpge_icpc_bulk_2vbb, "mul1/not_spms_coincident/is_ss_aoe/all")
```

A vector of simulations gives a dictionary of spectra, with the `simid`s as keys, ready
for [`build_forward_model`](@ref):

```julia
mc_dict = read_pdf(prod, [:hpge_icpc_bulk_2vbb, :liquid_argon_underground_Ar39];
                   groups = [:not_spms_coincident, :is_ss_aoe])
```
