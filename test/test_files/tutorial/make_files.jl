# Writes the fake input files of the tutorial in docs/src/tutorial.md:
# - one file for each simulated source with the histograms pdf/mul1/icpc and pdf/mul1/bege
#   and the number of simulated primaries,
# - data.lh5 with the measured histograms hist/mul1/icpc and hist/mul1/bege.
#
# Sources with a gamma line have a peak with a 1 keV width on a flat or linear
# continuum below the peak. The 2vbb spectrum is a broad Gaussian (600 keV width,
# centered at 1 MeV) and makes about 90% of the counts.
#
# Run with: julia --project=docs test/test_files/tutorial/make_files.jl
using LegendHDF5IO
using Distributions
using Random

Random.seed!(1)  # the same fake data in every build

dir = @__DIR__  # the files are written next to this script
lh5_binning = (binedges = (first = 0.0, last = 4000.0, step = 1.0), closedleft = true)
energy = collect(0.5:1.0:3999.5)  # bin centers in keV

# A peak with a 1 keV width containing 10% of the events, on a linear continuum
# 1 + slope * energy / peak up to the peak. The spectrum is normalised to 1.
function peak_shape(peak; slope = 0.0)
    continuum = (1 .+ slope .* energy ./ peak) .* (energy .< peak)
    line = exp.(-(energy .- peak) .^ 2 ./ 2)
    return 0.9 .* continuum ./ sum(continuum) .+ 0.1 .* line ./ sum(line)
end

# A broad Gaussian spectrum, normalised to 1
broad_shape(center, sigma) =
    (w = exp.(-(energy .- center) .^ 2 ./ (2 * sigma^2)); w ./ sum(w))

n_primaries = 1_000_000
efficiency = 0.3  # fraction of the primaries with an event in the detector
groups_lh5 = ["mul1/icpc", "mul1/bege"]

# name of the file => spectrum
sources = Dict(
    "cables_K40" => peak_shape(1461; slope = 0.0),
    "fibers_K40" => peak_shape(1461; slope = 1.0),
    "cables_Pb212_to_Pb208" => peak_shape(2615; slope = -0.5),
    "fibers_Pb212_to_Pb208" => peak_shape(2615; slope = 0.5),
    "hpge_icpc_bulk_2vbb" => broad_shape(1000, 600),
    "hpge_bege_bulk_2vbb" => broad_shape(950, 550),
)

for (name, spectrum) in sources
    lh5open(joinpath(dir, "$name.lh5"), "w") do f
        for g in groups_lh5
            f["pdf/$g/binning/axis_0"] = lh5_binning
            f["pdf/$g/weights"] = efficiency * n_primaries .* spectrum
        end
        f["nr_sim_events"] = n_primaries
    end
end

true_livetime = 1.0e6  # s
true_activities = (  # mBq
    cables_K40 = 200.0,
    fibers_K40 = 100.0,
    cables_Pb212_to_Pb208 = 50.0,
    fibers_Pb212_to_Pb208 = 30.0,
    two_vbb_hpge_bulk = 3500.0,
)
true_alpha = 2.0  # counts per 1 keV bin

two_vbb = 0.8 .* sources["hpge_icpc_bulk_2vbb"] .+ 0.2 .* sources["hpge_bege_bulk_2vbb"]

# the simulated spectra are summed over the two detector types, so the efficiency
# is counted twice
expected =
    1e-3 * true_livetime * 2 * efficiency .* (
        true_activities.cables_K40 .* sources["cables_K40"] .+
        true_activities.fibers_K40 .* sources["fibers_K40"] .+
        true_activities.cables_Pb212_to_Pb208 .* sources["cables_Pb212_to_Pb208"] .+
        true_activities.fibers_Pb212_to_Pb208 .* sources["fibers_Pb212_to_Pb208"] .+
        true_activities.two_vbb_hpge_bulk .* two_vbb
    ) .+ true_alpha

lh5open(joinpath(dir, "data.lh5"), "w") do f
    for g in groups_lh5
        f["hist/$g/binning/axis_0"] = lh5_binning
        f["hist/$g/weights"] = Float64.([rand(Poisson(l / 2)) for l in expected])
    end
end
