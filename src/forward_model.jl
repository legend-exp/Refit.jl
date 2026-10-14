"""
    build_forward_model(pdfs; e_low, e_high, blind, norm, binning)

Build a function that maps fit parameters to the predicted counts per energy bin.

`pdfs` is a dictionary `Dict(name => Histogram)` of simulated 1D energy spectra
that all share the same binning. Energies are in keV.
Only bins fully inside `[e_low, e_high]` and outside the `blind = (min, max)`
region are kept. `binning` merges that many consecutive bins first.

Returns `(forward_model, bin_edges)`:
- `bin_edges`: `(low, high)`, the lower and upper edges in keV of the kept bins
- `forward_model(params)`: predicted counts per kept bin for
  `params = (name1 = a1, name2 = a2, ...)`, with one activity per
  spectrum, in any order. Each spectrum is scaled by `activity * norm`. Use `norm`
  for unit conversions and the livetime. Each call allocates one array.

The time per call grows with the number of bins times the number of spectra.
"""
function build_forward_model(
    pdfs::AbstractDict{Symbol,<:Histogram{<:Real,1}};
    e_low::Real = 500,
    e_high::Real = 4000,
    blind::Union{Tuple{Real,Real},Nothing} = nothing,
    norm::Real = 1,
    binning::Int = 1,
)
    isempty(pdfs) && throw(ArgumentError("pdfs is empty"))
    names = Tuple(sort!(collect(keys(pdfs))))
    hists = map(k -> rebin_hist(pdfs[k], binning), names)

    edges = first(hists).edges[1]
    all(h -> h.edges[1] == edges, hists) || error("Histogram binning mismatch")

    is_valid = get_valid_bins(edges, (e_low, e_high), blind)
    bin_low = edges[1:(end-1)][is_valid]
    bin_high = edges[2:end][is_valid]
    mc = NamedTuple{names}(map(h -> h.weights[is_valid] .* norm, hists))
    # the same spectra as a vector, in the order of `names`. Julia unrolls a loop over
    # a tuple only up to about 30 entries; past that it falls back to a dynamic call per
    # spectrum, which costs more than adding the bins. A vector is looped over plainly.
    spectra = collect(values(mc))

    function forward_model(params::NamedTuple)
        # the activities in the order of the spectra. `keys(mc)` comes from the type of
        # `mc`, so the compiler sorts them here instead of on every call
        activities = NamedTuple{keys(mc)}(params)
        pred = zeros(length(bin_low))
        for j in eachindex(spectra)
            pred .+= getfield(activities, j) .* spectra[j]
        end
        return pred
    end
    return forward_model, (low = bin_low, high = bin_high)
end
