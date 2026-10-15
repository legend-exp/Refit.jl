# Writes the fake input files of the tutorial in docs/src/tutorial_fccd.md:
# - templates.lh5 with the generalised histogram pdf/geds, the spectra of a 60 keV
#   gamma line simulated on a grid of FCCD and transition layer (TL) thickness, and
#   the number of simulated primaries,
# - data.lh5 with the measured histogram hist/geds: the gamma line plus a flat
#   background.
#
# Each gamma deposits its 60 keV at a depth drawn uniformly from 0 to 2 mm below the
# surface of the detector. The charge collected from that depth is 0 in the dead
# layer (depth below FCCD - TL), rises linearly through the transition layer and is
# 1 below the FCCD. So the spectrum is a peak at 60 keV with probability
# (2 mm - FCCD) / 2 mm and a flat continuum from 0 to 60 keV with probability
# TL / 2 mm. Events from the dead layer have no energy and are not recorded. The
# measurement also has a flat background of 10 counts per keV from 0 to 100 keV.
#
# The histograms are written as structs with the fields of an LH5 histogram, and
# every group as a struct, so that the Python lh5 package can read the files too.
#
# Run with: julia --project=docs test/test_files/fccd/make_files.jl
using LegendHDF5IO
import LegendHDF5IO: Unitful
using Distributions
using StatsBase
using Random

dir = @__DIR__  # the files are written next to this script
edges = 0.0:1.0:100.0  # 1 keV bins

# fraction of the charge collected from the depth d, FCCD and TL in mm
function activeness(d, fccd, tl)
    d >= fccd && return 1.0
    d <= fccd - tl && return 0.0
    return (d - (fccd - tl)) / tl
end

# the counts of n gammas of 60 keV in 1 keV bins
function simulate(rng, n, fccd, tl)
    energy = 60 .* activeness.(2 .* rand(rng, n), fccd, tl)
    return fit(Histogram{Float64}, filter(>(0), energy), edges).weights
end

# an LH5 histogram with the bins of `edges`
lh5_histogram(counts) = (
    binning = (
        axis_0 = (
            binedges = (first = first(edges), last = last(edges), step = step(edges)),
            closedleft = true,
        ),
    ),
    weights = counts,
    isdensity = false,
)

# a struct of LH5 objects, from the names and values of its fields
lh5_struct(names, values) = NamedTuple{Tuple(Symbol.(names))}(Tuple(values))

# the templates, on a grid of 5 FCCD values times 6 TL values
rng = MersenneTwister(1)
n_primaries = 1_000_000
grid = (fccd = 0.6:0.2:1.4, tl = 0.0:0.1:0.5)  # mm

grid_info = lh5_struct(
    ["axis_$(d - 1)" for d = 1:length(grid)],
    [
        (
            name = String(name),
            values = (
                first = first(grid[name]) * Unitful.mm,
                last = last(grid[name]) * Unitful.mm,
                step = step(grid[name]) * Unitful.mm,
            ),
        ) for (d, name) in enumerate(keys(grid))
    ],
)
points = lh5_struct(
    ["fccd_$(i - 1)" for i in eachindex(grid.fccd)],
    [
        lh5_struct(
            ["tl_$(j - 1)" for j in eachindex(grid.tl)],
            [lh5_histogram(simulate(rng, n_primaries, fccd, tl)) for tl in grid.tl],
        ) for fccd in grid.fccd
    ],
)

lh5open(joinpath(dir, "templates.lh5"), "w") do f
    f["pdf"] = (geds = (grid_info = grid_info, grid = points),)
    f["nr_sim_events"] = n_primaries
end

# the measurement: 20000 gammas, a number known from the activity of the source and
# the livetime, and a flat background of 10 counts per keV
rng = MersenneTwister(2)
true_fccd, true_tl = 1.05, 0.27  # mm
data =
    simulate(rng, 20_000, true_fccd, true_tl) .+ rand(rng, Poisson(10.0), length(edges) - 1)

lh5open(joinpath(dir, "data.lh5"), "w") do f
    f["hist"] = (geds = lh5_histogram(Float64.(data)),)
end
