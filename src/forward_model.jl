"""
    build_forward_model(pdfs; e_low = 500, e_high = 4000, blind = nothing,
                        norm = 1, binning = 1)

Build a function that maps fit parameters to the predicted counts per energy bin.

`pdfs` is a dictionary `Dict(name => Histogram)` of simulated 1D energy spectra
that all share the same binning. Energies are in keV.
Only bins fully inside `[e_low, e_high]` and outside the `blind = (min, max)`
region are kept. `binning` merges that many consecutive bins first.

Returns `(forward_model, bin_edges)`:
- `bin_edges`: `(low, high)`, the lower and upper edges in keV of the kept bins
- `forward_model(params)`: predicted counts per kept bin for
  `params = (activities = (name1 = a1, name2 = a2, ...),)`, with one activity per
  spectrum, in any order. Each spectrum is scaled by `activity * norm`. Use `norm`
  for unit conversions and the livetime. Each call allocates one array.

The loop over spectra is unrolled at compile time, which only works for about 30
spectra or fewer. With more, each call makes extra small allocations and is slower.
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

    function forward_model(params::NamedTuple{(:activities,),<:Tuple{NamedTuple}})
        pred = zeros(length(bin_low))
        foreach(values(NamedTuple{keys(mc)}(params.activities)), values(mc)) do a, pdf
            pred .+= a .* pdf
        end
        return pred
    end
    return forward_model, (low = bin_low, high = bin_high)
end
