# Made-up files for the tutorial: simulated spectra and a measured spectrum in LH5 format.
using LegendHDF5IO
using Random
using Distributions

Random.seed!(1)

dir = mktempdir()
file_binning = (binedges = (first = 0.0, last = 4000.0, step = 1.0), closedleft = true)
energy = collect(0.5:1.0:3999.5)  # bin centers in keV

# Spectrum of one simulated source, normalised to 1: a peak with a 1 keV width, containing
# 10% of the events, on a linear continuum up to the peak. The continuum is
# 1 + slope * energy / peak, so a slope of 0 is flat.
function peak_shape(peak; slope = 0.0)
    continuum = (1 .+ slope .* energy ./ peak) .* (energy .< peak)
    continuum ./= sum(continuum)
    line = exp.(-(energy .- peak) .^ 2 ./ (2 * 1.0^2))
    return 0.9 .* continuum .+ 0.1 .* line ./ sum(line)
end

# A broad Gaussian spectrum with no peak, normalised to 1.
function broad_shape(center, sigma)
    w = exp.(-(energy .- center) .^ 2 ./ (2 * sigma^2))
    return w ./ sum(w)
end

n_primaries = 1_000_000
efficiency = 0.3  # fraction of the primaries with an event in the detector

# name of the pdf file => spectrum
sources = Dict(
    "cables_K40" => peak_shape(1461; slope = 0.0),
    "fibers_K40" => peak_shape(1461; slope = 1.0),
    "cables_Pb212_to_Pb208" => peak_shape(2615; slope = -0.5),
    "fibers_Pb212_to_Pb208" => peak_shape(2615; slope = 0.5),
    "hpge_icpc_bulk_2vbb" => broad_shape(1000, 600),
    "hpge_bege_bulk_2vbb" => broad_shape(950, 550),
)

groups = ["mul1/icpc", "mul1/bege"]  # the two detector types of the spectrum

for (name, spectrum) in sources
    lh5open(joinpath(dir, "$name.lh5"), "w") do f
        for g in groups
            f["pdf/$g/binning/axis_0"] = file_binning
            f["pdf/$g/weights"] = efficiency * n_primaries .* spectrum
        end
        f["nr_sim_events"] = n_primaries
    end
end

# true activities in mBq, and the livetime in s
true_activities = (
    cables_K40 = 200.0,
    fibers_K40 = 100.0,
    cables_Pb212_to_Pb208 = 50.0,
    fibers_Pb212_to_Pb208 = 30.0,
    two_vbb_hpge_bulk = 8.0,
)
true_alpha = 2.0  # counts per 1 keV bin
true_livetime = 1.0e6

# The measured spectrum has Poisson fluctuations. The simulated spectra are the sum
# of the two detector types, so the efficiency is counted twice.
two_vbb = 0.8 .* sources["hpge_icpc_bulk_2vbb"] .+ 0.2 .* sources["hpge_bege_bulk_2vbb"]
expected =
    1e-3 * true_livetime * 2 * efficiency .* (
        true_activities.cables_K40 .* sources["cables_K40"] .+
        true_activities.fibers_K40 .* sources["fibers_K40"] .+
        true_activities.cables_Pb212_to_Pb208 .* sources["cables_Pb212_to_Pb208"] .+
        true_activities.fibers_Pb212_to_Pb208 .* sources["fibers_Pb212_to_Pb208"] .+
        true_activities.two_vbb_hpge_bulk .* two_vbb
    ) .+ true_alpha

lh5open(joinpath(dir, "data.lh5"), "w") do f
    for g in groups
        f["hist/$g/binning/axis_0"] = file_binning
        f["hist/$g/weights"] = Float64.([rand(Poisson(l / 2)) for l in expected])
    end
end
