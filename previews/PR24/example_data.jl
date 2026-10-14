# Made-up spectra used by the examples in the user manual.
using Refit
using StatsBase

edges = 0.0:10.0:4000.0  # 10 keV bins
energy = collect(edges[1:(end-1)]) .+ 5.0  # bin centers in keV

# One simulated spectrum: counts per decay in each bin, with an efficiency of 0.5
function simulated_spectrum(continuum, peaks)
    peak(e0, h) = h .* exp.(-(energy .- e0) .^ 2 ./ (2 * 6.0^2))  # 6 keV width
    shape =
        continuum .+ sum((peak(e0, h) for (e0, h) in peaks); init = zeros(length(energy)))
    return Histogram(edges, 0.5 .* shape ./ sum(shape))
end

mc_dict = Dict(
    :K40_hpge_bulk => simulated_spectrum(exp.(-energy ./ 800), [(1461, 3.0)]),
    :Pb212_to_Pb208_hpge_bulk =>
        simulated_spectrum(exp.(-energy ./ 1500), [(583, 2.0), (2615, 4.0)]),
    :two_vbb_hpge_bulk => simulated_spectrum(energy .* exp.(-energy ./ 500), []),
)

# Measured spectrum: made-up activities in mBq and a livetime of 3 million seconds
livetime = 3.0e6
true_activities =
    (K40_hpge_bulk = 40.0, Pb212_to_Pb208_hpge_bulk = 10.0, two_vbb_hpge_bulk = 100.0)
data = Histogram(
    edges,
    round.(
        sum(
            1e-3 * livetime * true_activities[k] .* mc_dict[k].weights for
            k in keys(true_activities)
        ),
    ),
)
