# Made-up files for the tutorial: simulated spectra and a measured spectrum in LH5 format.
using LegendHDF5IO
using Random
using Distributions

Random.seed!(1)

dir = mktempdir()
file_binning = (binedges = (first = 0.0, last = 4000.0, step = 1.0), closedleft = true)
energy = collect(0.5:1.0:3999.5)  # bin centers in keV

# Counts per bin of one simulated source with a smooth continuum and an optional peak.
# The spectrum is normalised to 1, so that its sum is the fraction of the primaries
# that give an event.
function shape(slope; peak = nothing)
    w = exp.(-energy ./ slope)
    w ./= sum(w)
    if !isnothing(peak)
        p = exp.(-(energy .- peak) .^ 2 ./ (2 * 2.0^2))  # 2 keV width
        w = 0.9 .* w .+ 0.1 .* p ./ sum(p)
    end
    return w
end

n_primaries = 1_000_000
efficiency = 0.3  # fraction of the primaries with an event in the detector

# name of the pdf file => spectrum
sources = Dict(
    "cables_K40" => shape(900; peak = 1461),
    "fibers_K40" => shape(700; peak = 1461),
    "cables_Pb212_to_Pb208" => shape(1200; peak = 2615),
    "fibers_Pb212_to_Pb208" => shape(1000; peak = 2615),
    "hpge_icpc_bulk_2vbb" => shape(500),
    "hpge_bege_bulk_2vbb" => shape(450),
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
