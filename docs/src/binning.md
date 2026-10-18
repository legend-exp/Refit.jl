# Binning
`Refit.jl` contains several options for the binning used in the binned likelihood.
The keyword `binning` of [`build_forward_model`](@ref) rebins the histograms before the
energy range is applied.

```@setup demo
include("example_data.jl")
```

## Types of binning

There are three options:

| `binning` | Bins |
|---|---|
| `20` | every 20 consecutive bins merged into one |
| `[500, 520, 600, 1000]` or `500:20:3000` | these bin edges in keV, with bins of any width |
| `[500:20:1400, 1500:50:3000]` | several ranges, each with its own edges |

With several ranges, the regions between them are skipped, here 1400 to 1500 keV. Ranges
that touch are joined. Every edge must be an edge of the histograms, and each new bin is
the sum of the bins inside it. A range can also be a list of edges, as the last one in
this example, which skips 1400 to 1500 keV and 2500 to 2600 keV:

```@example demo
model_r, bin_edges_r = build_forward_model(
    mc_dict;
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,
    binning = [500:20:1400, 1500:100:2500, [2600, 2610, 2620, 2630]],
)

# the bins around the first skipped region, lower and upper edges in keV
bin_edges_r.low[44:47], bin_edges_r.high[44:47]
```

With bins of different widths, the counts per bin are not comparable from one bin to the
next. Divide them by the width of the bin to compare or plot them as counts per keV.

## Blinded regions

A region is blinded by leaving it out between two ranges. Here the bins are 10 keV wide
from 500 to 3200 keV, and 2010 to 2070 keV, around the Q-value of ``^{76}``Ge at
2039 keV, is not fitted:

```@example demo
# two ranges, with the blinded region between them
blinded = [500:10:2010, 2070:10:3200]
model_b, bin_edges_b = build_forward_model(
    mc_dict;
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,
    binning = blinded,
)

# the last bin below and the first bin above the blinded region, in keV
i = findfirst(==(2010.0), bin_edges_b.high)
(bin_edges_b.low[i], bin_edges_b.high[i]), (bin_edges_b.low[i+1], bin_edges_b.high[i+1])
```

[`get_priors`](@ref) takes the same `binning`, so it does not look at the blinded data
either:

```@example demo
priors_b = get_priors(
    mc_dict,
    data;
    e_low = 500,
    e_high = 3200,
    norm = 1e-3 * livetime,
    binning = blinded,
)
```

## Binning from a YAML file

In a YAML file a range is written as a string `"start:step:stop"`:

```yaml
binning:
  - 500:20:1400
  - 1500:100:2500
  - [2600, 2610, 2620, 2630]
```

```julia
using YAML

binning = YAML.load_file("fit.yaml")["binning"]
model, bin_edges = build_forward_model(mc_dict; binning, e_low = 500, e_high = 3200)
```

On a single line, the ranges need quotes: `["500:20:1400", "1500:100:2500"]`.

!!! warning "The end of a range"
    A string always ends at its `stop`. If the steps do not reach it, the last bin is
    narrower: `"500:200:2015"` has the edges 500, 700, ..., 1900 and 2015 keV. A Julia
    range follows the rule of Julia instead and stops at its last step: `500:200:2015`
    is `500:200:1900`, and the region from 1900 to 2015 keV is then not in the range, so
    it is skipped without a warning.

## Uniform binning with peaks

Refit.jl provides a convenience function to help build a standard type of binning.
[`uniform_binning_with_peaks`](@ref) gives uniform bins, except that each peak region is
a single bin. The bins next to a peak are between one and two bin widths wide, and every
edge is an edge of the histogram given:

```@example demo
# 50 keV bins, with the lines of K-40 and Tl-208 each as one bin
peaks = [(1450, 1470), (2600, 2630)]  # keV, on edges of the 10 keV bins
binning_p = uniform_binning_with_peaks(mc_dict[:K40_hpge_bulk], 50; peaks, low = 500, high = 3200)
binning_p[18:22]
```

With `blinded = (min, max)` the region is left out instead of becoming a bin. The result
is two ranges of edges, below and above it:

```@example demo
binning_b = uniform_binning_with_peaks(
    mc_dict[:K40_hpge_bulk], 50; peaks, blinded = (2010, 2070), low = 500, high = 3200,
)
last(binning_b[1]), first(binning_b[2])
```
