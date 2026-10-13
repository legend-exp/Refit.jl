# Likelihood

`poisson_likelihood` compares measured and predicted counts in each bin. The counts
in a bin follow a Poisson distribution with the predicted counts ``\lambda_i`` as
mean, so for measured counts ``k_i`` the log-likelihood is

```math
\ln L = \sum_i \left( k_i \ln \lambda_i - \lambda_i \right).
```

The constant term ``-\ln k_i!`` is left out. It does not depend on the parameters,
so it does not change the position of the maximum or the shape of the posterior.

## Measured counts

The likelihood compares two flat vectors of the same length: the measured counts and
the prediction of the [forward model](forward_model.md). Entry `i` of both is the
number of counts in the same energy bin. The measured counts are **not** a histogram.

`get_counts` extracts them from the data histogram `data`. It takes the lower and
upper edge of each bin, in keV:

```@setup demo
include("example_data.jl")
```

```@example demo
# 100 keV bins between 500 and 3200 keV
low = collect(500.0:100.0:3100.0)  # lower edge of each bin
bin_edges = (low = low, high = low .+ 100.0);  # upper edge = lower edge + 100 keV

# measured counts in each bin, as a Vector{Float64}
counts = get_counts(data, bin_edges)
```

The bins can be any bins whose edges are also edges of `data`, otherwise an error is
raised. A bin can span several bins of `data`, which are then summed.

The bin edges can also come from the forward model. `build_forward_model` returns the
edges of the bins it predicts, whatever energy range, blinding and binning it uses.
Using them guarantees that the measured and predicted bins match:

```@example demo
# merge pairs of bins (20 keV bins) and keep 500 to 3200 keV
model, bin_edges = build_forward_model(
    mc_dict;
    e_low = 500,
    e_high = 3200,
    binning = 2,
    norm = 1e-3 * livetime,
);

# the measured counts in exactly the bins of the model
counts = get_counts(data, bin_edges)
counts
```

## Building the likelihood

```@example demo
using DensityInterface

# the likelihood of the measured counts as a function of the parameters
likelihood = build_likelihood(counts, model);

params = (
    K40_hpge_bulk = 40.0,
    Pb212_to_Pb208_hpge_bulk = 10.0,
    two_vbb_hpge_bulk = 100.0,
)

# log-likelihood for these parameters
logdensityof(likelihood, params)
```

The value is largest when the prediction matches the data. Here the data were made
with these activities, so changing one of them lowers the value:

```@example demo
other = merge(params, (K40_hpge_bulk = 50.0,))
logdensityof(likelihood, other)
```

`build_likelihood` handles one measured spectrum and one model. To fit several
spectra together, add the log-likelihoods of each one.

`poisson_likelihood(counts, pred)` is the underlying function. It raises a
`DimensionMismatch` if the two vectors differ in length. A bin with a predicted
count of zero is allowed: a tiny number (`eps()`) is added to the prediction inside
the logarithm, so the result stays finite.
