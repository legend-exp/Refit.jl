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
    read_pdf(file; groups = ["pdf/mul_surv"], n_prim_field = "nr_sim_events")

Read a simulated spectrum from `file` as counts per simulated primary.

The histograms in `groups` are summed, then divided by the number of simulated primaries
stored in the field `n_prim_field`. To scale by anything else, change the weights of the
result: `h.weights ./= factor`.
"""
function read_pdf(
    file::String;
    groups::Vector{String} = ["pdf/mul_surv"],
    n_prim_field::String = "nr_sim_events",
)
    h = add_histograms([read_hist(g, file) for g in groups])
    h.weights ./= read_scalar(n_prim_field, file)
    return h
end
