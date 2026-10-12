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
