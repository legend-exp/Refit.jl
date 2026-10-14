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

The spectra are stored in one table, one column per spectrum, and added in a plain
loop. The time per call grows with the number of bins times the number of spectra.
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
    # one column per spectrum, already multiplied by norm
    table = zeros(length(bin_low), length(names))
    for (j, h) in enumerate(hists)
        table[:, j] .= h.weights[is_valid] .* norm
    end
    column_of = Dict{Symbol,Int}(name => j for (j, name) in enumerate(names))

    forward_model(params::NamedTuple) =
        add_spectra!(zeros(length(bin_low)), table, column_of, params)
    return forward_model, (low = bin_low, high = bin_high)
end

# Adds `params[name] * column` for each column of the table. Looking up a parameter by
# name is slow when there are many, so the parameters are visited by position instead.
function add_spectra!(pred, table, column_of::Dict{Symbol,Int}, params::NamedTuple)
    found = 0
    for (pos, name) in enumerate(keys(params))
        j = get(column_of, name, 0)
        j == 0 && continue
        found += 1
        a = getfield(params, pos)
        @inbounds @simd for i in eachindex(pred)
            pred[i] += a * table[i, j]
        end
    end
    if found < length(column_of)  # raises the error for the missing parameter
        foreach(name -> params[name], keys(column_of))
    end
    return pred
end
