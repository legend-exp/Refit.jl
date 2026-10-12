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

`get_counts` extracts them from the data histogram, using the `bin_edges` returned
by `build_forward_model`. This guarantees that the bins match, whatever energy range,
blinding and binning the model uses.

```julia
model, bin_edges = build_forward_model(mc_dict; e_low = 500, e_high = 3200, binning = 2)

counts = get_counts(data, bin_edges)  # Vector{Float64}, one entry per bin of the model
```

`data` is a `Histogram` of the measured spectrum. The edges in `bin_edges` must be
edges of `data`, otherwise an error is raised. If the model merges bins with
`binning`, the data bins are summed in the same way.

## Building the likelihood

```julia
likelihood = build_likelihood(counts, model)

using DensityInterface
logdensityof(likelihood, params)  # params as for the model
```

`build_likelihood` handles one measured spectrum and one model. To fit several
spectra together, add the log-likelihoods of each one.

`poisson_likelihood(counts, pred)` is the underlying function. It raises a
`DimensionMismatch` if the two vectors differ in length. A bin with a predicted
count of zero is allowed: a tiny number (`eps()`) is added to the prediction inside
the logarithm, so the result stays finite.
