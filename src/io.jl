using LegendHDF5IO

"""
     read_hist(group, file)

Read a histogram in group `group` from an lh5 file `file` and return a `Histogram` object.

The histogram format is described in the LH5 specification. The `group` should contain a `binning` dataset and a `weights` dataset.
"""
function read_hist(group::String, file::String)

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
