using LegendHDF5IO
using YAML

"""
     read_hist(group, file)

Read a histogram in group `group` from an lh5 file `file` and return a `Histogram` object.

The histogram format is described in the LH5 specification. The `group` should contain a `binning` dataset and a `weights` dataset. Histograms written from Julia (for example by [`write_best_fit`](@ref)) are read too.
"""
function read_hist(group::String, file::String)

    # files written from Julia count the axes from 1, files from the simulation from 0
    from_julia = lh5open(f -> haskey(f, "$group/binning/axis_1"), "$file", "r")
    if from_julia
        # read the edges and weights into memory before the file is closed
        return lh5open("$file", "r") do f
            h = f[group]
            Histogram(map(collect, h.edges), collect(h.weights))
        end
    end

    binning = lh5open("$file", "r") do f
        f["$group/binning/axis_0"]
    end
    weights = lh5open("$file", "r") do f
        f["$group/weights"][:]
    end
    h = LegendHDF5IO._nt_to_histogram((
        binning = [binning],
        weights = weights,
        isdensity = false,
    ))

    return h

end


"""
    read_scalar(field, file)

Read a Scalar from `file` called `field`
"""
function read_scalar(field::String, file::String)
    lh5open("$file", "r") do f
        f[field]
    end
end

"""
    read_pdf(file; groups = "pdf/mul_surv", weights = nothing, n_prim_field = "nr_sim_events")

Read a simulated spectrum from `file`, by default as counts per simulated primary.

`groups` is the name of one group, or a vector of names. The histograms in these groups
are multiplied by `weights` (a vector with one factor per group, by default 1) and summed,
then divided by the number of simulated primaries stored in the field `n_prim_field`.
With `n_prim_field = nothing` the sum is not divided.
To scale by anything else, change the weights of the result: `h.weights ./= factor`.
"""
function read_pdf(
    file::String;
    groups::Union{String,AbstractVector{String}} = "pdf/mul_surv",
    weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    n_prim_field::Union{String,Nothing} = "nr_sim_events",
)
    groups = vcat(groups)
    weights = isnothing(weights) ? ones(length(groups)) : weights
    length(weights) == length(groups) ||
        throw(DimensionMismatch("groups and weights differ in length"))

    hists = [read_hist(g, file) for g in groups]
    for (h, w) in zip(hists, weights)
        h.weights .*= w
    end
    total = add_histograms(hists)
    isnothing(n_prim_field) || (total.weights ./= read_scalar(n_prim_field, file))
    return total
end

"""
    read_pdf(files::AbstractDict; groups = "pdf/mul_surv", n_prim_field = "nr_sim_events")

Read the same source simulated in several files, for example once per detector type,
and return the weighted average of the spectra.

`files` maps each file name to its weight, for example the mass of the detector type.
The weights are divided by their sum. Each file is read with the single-file method,
so it is divided by its own number of primaries.
"""
function read_pdf(
    files::AbstractDict{String,<:Real};
    groups::Union{String,AbstractVector{String}} = "pdf/mul_surv",
    n_prim_field::Union{String,Nothing} = "nr_sim_events",
)
    total = sum(values(files))
    hists = map(collect(files)) do (file, w)
        h = read_pdf(file; groups, n_prim_field)
        h.weights .*= w / total
        h
    end
    return add_histograms(hists)
end

"""
    read_generalised_hist(group, file)

Read the [`GeneralisedHistogram`](@ref) in `group` of the LH5 file `file`.

The group holds the description of the grid in `grid_info` and the histogram of every
grid point in `grid`, for example for a grid of FCCD and transition layer thickness:

```text
group
├── grid_info
│   ├── axis_0
│   │   ├── name      "fccd"
│   │   └── values    first, last, step (with the units as an attribute)
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

Each `axis_d` is one parameter of the grid: its name and evenly spaced values from
`first` to `last` in steps of `step`. The groups in `grid` are nested in the order of
the axes, and `<name>_<i>` holds the `i`-th value of that parameter, counted from 0.
Every combination must be present, with the same binning. The values are read in the
units of the file, which are not kept.
"""
function read_generalised_hist(group::String, file::String)
    names, values = lh5open(file, "r") do f
        n_axes = length(keys(f["$group/grid_info"]))
        names = Symbol[]
        values = AbstractRange[]
        for d = 0:(n_axes-1)
            axis = "$group/grid_info/axis_$d"
            haskey(f, axis) || throw(ArgumentError("$file has no $axis"))
            push!(names, Symbol(f["$axis/name"]))
            # strip the units, if any: a quantity divided by its unit is a number
            first, last, step = (
                (x = f["$axis/values/$k"]; x / oneunit(x)) for
                k in ("first", "last", "step")
            )
            n = round(Int, (last - first) / step) + 1
            r = range(first; step, length = n)
            Base.last(r) ≈ last || throw(
                ArgumentError("$axis does not go from $first to $last in steps of $step"),
            )
            push!(values, r)
        end
        return names, values
    end

    # the path of the histogram at each grid point, for example "fccd_2/tl_0"
    paths = map(CartesianIndices(Tuple(length.(values)))) do idx
        join(("$(names[d])_$(idx[d] - 1)" for d in eachindex(names)), "/")
    end
    lh5open(file, "r") do f
        for p in paths
            haskey(f, "$group/grid/$p") ||
                throw(ArgumentError("$file has no $group/grid/$p"))
        end
    end
    hists = map(p -> read_hist("$group/grid/$p", file), paths)
    return GeneralisedHistogram(hists; NamedTuple{Tuple(names)}(Tuple(values))...)
end

"""
    read_generalised_pdf(file; groups = "pdf/mul_surv", weights = nothing, n_prim_field = "nr_sim_events")

Read a simulated spectrum on a grid of parameters from `file`, by default as counts per
simulated primary. The same as [`read_pdf`](@ref), for a [`GeneralisedHistogram`](@ref):
the generalised histograms in `groups` are multiplied by `weights` and summed, grid point
by grid point, then divided by the number of simulated primaries in the field
`n_prim_field`. All grid points must have the same number of simulated primaries.
"""
function read_generalised_pdf(
    file::String;
    groups::Union{String,AbstractVector{String}} = "pdf/mul_surv",
    weights::Union{Nothing,AbstractVector{<:Real}} = nothing,
    n_prim_field::Union{String,Nothing} = "nr_sim_events",
)
    groups = vcat(groups)
    weights = isnothing(weights) ? ones(length(groups)) : weights
    length(weights) == length(groups) ||
        throw(DimensionMismatch("groups and weights differ in length"))

    ghs = [read_generalised_hist(g, file) for g in groups]
    for (gh, w) in zip(ghs, weights)
        foreach(c -> c .*= w, gh.counts)
    end
    total = add_histograms(ghs)
    if !isnothing(n_prim_field)
        n_primaries = read_scalar(n_prim_field, file)
        foreach(c -> c ./= n_primaries, total.counts)
    end
    return total
end

"""
    write_best_fit(file, hists; data = nothing)

Write the histograms of [`get_best_fit_histograms`](@ref) to the lh5 file `file`: the
total in `best_fit/total` and each source in `best_fit/<name>`. An existing file is
overwritten.

With `data` (a `Histogram`, for example from [`counts_to_histogram`](@ref)) the
measured spectrum is also written, in `data`. The histograms can be read back with
[`read_hist`](@ref).
"""
function write_best_fit(
    file::AbstractString,
    hists::AbstractDict{Symbol,<:Histogram{<:Real,1}};
    data::Union{Nothing,Histogram{<:Real,1}} = nothing,
)
    lh5open(file, "w") do f
        for (name, h) in hists
            f["best_fit/$name"] = h
        end
        isnothing(data) || (f["data"] = data)
    end
    return file
end

"""
    write_parameter_values(file, parameters, intervals; level = 0.68)

Write the values of the parameters, for example the activities, to the YAML file `file`.
For each parameter it saves the best fit value from `parameters` and the interval `low` to `high` from `intervals`. Both are
named tuples with the same names, and each interval is an `Interval`.

The intervals should be the smallest credible intervals, the shortest intervals that
contain the probability `level` of the posterior. BAT.jl computes them:
`map(first, BAT.smallest_credible_intervals(samples; p = 0.68, mode = :connected))`.
The `level` is saved in the file as `credible_level`.
"""
function write_parameter_values(
    file::AbstractString,
    parameters::NamedTuple,
    intervals::NamedTuple;
    level::Real = 0.68,
)
    # ordered dictionaries keep the order of the parameters in the file
    values = OrderedDict(
        string(name) => OrderedDict(
            "best_fit" => Float64(parameters[name]),
            "low" => Float64(intervals[name].left),
            "high" => Float64(intervals[name].right),
        ) for name in keys(parameters)
    )
    YAML.write_file(
        file,
        OrderedDict("credible_level" => Float64(level), "parameters" => values),
    )
    return file
end
