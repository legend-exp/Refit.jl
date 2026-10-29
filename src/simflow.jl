"""
    SimProd(path)

A production cycle of [legend-simflow](https://github.com/legend-exp/legend-simflow) in
the directory `path`, the one with the `simflow-config.yaml` file.

The experiment name and the directory of the `pdf` tier files are read from
`simflow-config.yaml`. The `pdf` tier has one file per simulation (`simid`), named
`<experiment>-<simid>-tier_pdf.lh5`.

The list of the simulations and the list of the histograms are read from disk the first
time they are needed, and kept for later calls. The histograms are read from the file of
the first simulation: all files of a cycle must have the same histograms. Files added
to the cycle later are seen by a new `SimProd`.

# Examples
```julia
prod = SimProd("prod/l1000dsg01/v1.1.0")
simids(prod)  # the simulations in the cycle
pdf_file(prod, :hpge_icpc_bulk_2vbb)  # the path of the pdf tier file
```
"""
mutable struct SimProd
    const path::String
    const experiment::String
    const pdf_dir::String
    # filled at the first use, see _simids and _pdf_groups
    simids::Union{Nothing,Vector{Symbol}}
    pdf_groups::Union{Nothing,Vector{String}}
end

function SimProd(path::AbstractString)
    path = abspath(path)
    config_file = joinpath(path, "simflow-config.yaml")
    isfile(config_file) || throw(ArgumentError("$path has no simflow-config.yaml"))
    config = YAML.load_file(config_file)
    # "$_" in the paths of the configuration stands for the cycle directory
    pdf_dir = normpath(replace(config["paths"]["tier"]["pdf"], raw"$_" => path))
    return SimProd(path, config["experiment"], pdf_dir, nothing, nothing)
end

"""The simids of `prod`, read from disk at the first call."""
function _simids(prod::SimProd)
    if isnothing(prod.simids)
        prefix, suffix = "$(prod.experiment)-", "-tier_pdf.lh5"
        files =
            filter(f -> startswith(f, prefix) && endswith(f, suffix), readdir(prod.pdf_dir))
        prod.simids = [Symbol(chopsuffix(chopprefix(f, prefix), suffix)) for f in files]
    end
    return prod.simids
end

"""Paths of the histograms below the HDF5 group `group`, relative to it."""
function _hist_paths(group::HDF5.Group)
    haskey(group, "binning") && return [""]
    paths = String[]
    for k in keys(group)
        child = group[k]
        child isa HDF5.Group || continue
        append!(paths, (isempty(p) ? k : "$k/$p" for p in _hist_paths(child)))
    end
    return paths
end

"""The histogram paths of `prod`, read from the file of the first simulation at the first call."""
function _pdf_groups(prod::SimProd)
    if isnothing(prod.pdf_groups)
        file = pdf_file(prod, first(_simids(prod)))
        prod.pdf_groups = HDF5.h5open(f -> _hist_paths(f["pdf"]), file)
    end
    return prod.pdf_groups
end

"""Each histogram path of `prod` with its categories and HPGe group, if it has both."""
function _cuts_and_hpge(prod::SimProd)
    parts = (p => Symbol.(split(p, '/')) for p in _pdf_groups(prod))
    return [
        (path = p, cuts = ps[1:(end-1)], hpge = ps[end]) for
        (p, ps) in parts if length(ps) > 1
    ]
end

"""
    simids(prod::SimProd)

The simulations (`simid`) of `prod` with a `pdf` tier file, as `Symbol`s in alphabetical
order.
"""
simids(prod::SimProd) = copy(_simids(prod))

"""
    pdf_file(prod::SimProd, simid::Symbol)

The path of the `pdf` tier file of the simulation `simid` in `prod`.
"""
function pdf_file(prod::SimProd, simid::Symbol)
    simid in _simids(prod) ||
        throw(ArgumentError("$(prod.path) has no pdf tier file for $simid"))
    return joinpath(prod.pdf_dir, "$(prod.experiment)-$simid-tier_pdf.lh5")
end

"""
    pdf_groups(prod::SimProd)

The paths of all histograms below `pdf/` in the `pdf` tier of `prod`, for example
`"mul1/not_spms_coincident/is_ss_aoe/all"`. Only the file of the first simulation is
read, the first time: all files of a cycle have the same histograms.
"""
pdf_groups(prod::SimProd) = copy(_pdf_groups(prod))

"""
    cut_groups(prod::SimProd)

The combinations of event categories in the `pdf` tier of `prod`, one for each
histogram of a HPGe group, for example `[:is_ss_aoe, :mul1, :not_spms_coincident]`.
The categories of a combination are in alphabetical order: their order has no meaning.
[`read_pdf`](@ref) sums the combinations that contain all the categories of `groups`.
"""
cut_groups(prod::SimProd) = unique(sort(h.cuts) for h in _cuts_and_hpge(prod))

"""
    hpge_groups(prod::SimProd)

The HPGe groups of the histograms in the `pdf` tier of `prod`, for example `:all`.
Any of them can be given as `hpge_group` to [`read_pdf`](@ref).
"""
hpge_groups(prod::SimProd) = unique(h.hpge for h in _cuts_and_hpge(prod))

Base.show(io::IO, prod::SimProd) = print(io, "SimProd(", repr(prod.path), ")")

function Base.show(io::IO, ::MIME"text/plain", prod::SimProd)
    ids = _simids(prod)
    println(io, "SimProd at ", prod.path)
    println(io, "├─ experiment: ", prod.experiment)
    isempty(ids) && return print(io, "└─ 0 simulations")
    println(io, "├─ $(length(ids)) simulations:")
    # in the REPL, only the first and last few of a long list, like Julia's arrays
    names = string.(ids)
    if get(io, :limit, false) && length(names) > 5
        names = [names[1:3]; "⋮"; names[end]]
    end
    foreach(n -> println(io, "│    ", n), names)
    println(io, "├─ cut groups:")
    foreach(c -> println(io, "│    ", join(c, " ∧ ")), cut_groups(prod))
    print(io, "└─ HPGe groups: ", join(hpge_groups(prod), ", "))
end

"""The entries of `_cuts_and_hpge` that match `groups` and `hpge_group`, `nothing` matching any group."""
function _matching_hists(prod, groups, hpge_group)
    # one list of cuts that must all apply, or a vector of such lists of which any applies
    lists =
        groups isa Symbol ? [[groups]] : all(g -> g isa Symbol, groups) ? [groups] : groups
    hpge = isnothing(hpge_group) ? nothing : vcat(hpge_group)
    return filter(_cuts_and_hpge(prod)) do h
        (isnothing(hpge) || h.hpge in hpge) && any(list -> list ⊆ h.cuts, lists)
    end
end

"""
    has_pdf(prod::SimProd, simid::Symbol; groups = Symbol[], hpge_group = nothing)
    has_pdf(prod::SimProd, simid::Symbol, path)

Whether `prod` has the spectrum of the simulation `simid`, without reading it. After the
first call, it does not read the disk (see [`SimProd`](@ref)).

With only `simid`, whether the simulation has a `pdf` tier file. With `groups` or
`hpge_group`, whether [`read_pdf`](@ref) finds a histogram for them: they have the same
meaning as there. An empty `groups` has no cuts and matches any histogram, and
`hpge_group = nothing` matches any HPGe group. With a vector of HPGe groups, whether
each of them has a histogram, so an empty vector is always `true`. With `path`, a
string, whether the file has a histogram at this path below `pdf/`, as listed by
[`pdf_groups`](@ref).

# Examples
```julia
has_pdf(prod, :hpge_icpc_bulk_2vbb)  # the simulation is in the cycle
has_pdf(prod, :hpge_icpc_bulk_2vbb; groups = [:is_ss_aoe, :not_spms_coincident], hpge_group = :all)
has_pdf(prod, :hpge_icpc_bulk_2vbb; hpge_group = :icpc)  # in any category
has_pdf(prod, :hpge_icpc_bulk_2vbb, "mul1/not_spms_coincident/is_ss_aoe/all")
```
"""
function has_pdf(
    prod::SimProd,
    simid::Symbol;
    groups::Union{Symbol,AbstractVector} = Symbol[],
    hpge_group::Union{Nothing,Symbol,AbstractVector} = nothing,
)
    simid in _simids(prod) || return false
    # without a selection, do not read the histogram paths from disk
    groups isa AbstractVector && isempty(groups) && isnothing(hpge_group) && return true
    found = [h.hpge for h in _matching_hists(prod, groups, hpge_group)]
    return isnothing(hpge_group) ? !isempty(found) : vcat(hpge_group) ⊆ found
end

function has_pdf(prod::SimProd, simid::Symbol, path::AbstractString)
    return simid in _simids(prod) && strip(path, '/') in _pdf_groups(prod)
end

"""
    read_pdf(prod::SimProd, simid::Symbol; groups = :mul1, hpge_group = :all)
    read_pdf(prod::SimProd, simids::AbstractVector{Symbol}; groups = :mul1, hpge_group = :all)

Read the spectrum of the simulation `simid` in `prod`, as counts per simulated primary.
With a vector of `simids`, return a dictionary from each `simid` to its spectrum, which
can be given to [`build_forward_model`](@ref).

The `pdf` tier sorts the events into categories: by multiplicity, the LAr veto and the
A/E cut, for example. It stores one histogram for each combination of categories (see
[`cut_groups`](@ref)) and each HPGe group. These histograms do not overlap, so any
event selection is a sum of them.

`groups` is a category or a vector of categories, in any order. It selects the
histograms whose combination contains all these categories. A vector of such vectors
selects the histograms matched by any of them, each counted once. `hpge_group` is the
HPGe group, or a vector of HPGe groups to sum the histograms of all of them.

# Examples
```julia
read_pdf(prod, :hpge_icpc_bulk_2vbb)  # all multiplicity-1 events
read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = :not_spms_coincident)  # after the LAr veto
read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = :is_ss_aoe)  # after the A/E cuts
# after the LAr veto and the A/E cuts, in any order
read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = [:is_ss_aoe, :not_spms_coincident])
# after the low-side A/E cut only
read_pdf(prod, :hpge_icpc_bulk_2vbb; groups = [[:is_ss_aoe], [:is_high_aoe]])
```
"""
function read_pdf(
    prod::SimProd,
    simid::Symbol;
    groups::Union{Symbol,AbstractVector} = :mul1,
    hpge_group::Union{Symbol,AbstractVector} = :all,
)
    file = pdf_file(prod, simid)
    paths = [h.path for h in _matching_hists(prod, groups, hpge_group)]
    isempty(paths) && throw(
        ArgumentError("no histogram of HPGe group $hpge_group matches $groups in $file"),
    )
    return read_pdf(file; groups = ["pdf/$p" for p in paths])
end

function read_pdf(prod::SimProd, simids::AbstractVector{Symbol}; kwargs...)
    return OrderedDict(s => read_pdf(prod, s; kwargs...) for s in simids)
end
