# Likelihood

`poisson_likelihood` compares measured and predicted counts in each bin. The counts
in a bin follow a Poisson distribution with the predicted counts ``\lambda_i`` as
mean, so for measured counts ``k_i`` the log-likelihood is

```math
\ln L = \sum_i \left( k_i \ln \lambda_i - \lambda_i \right).
```

The constant term ``-\ln k_i!`` is left out. It does not depend on the parameters,
so it does not change the position of the maximum or the shape of the posterior.

## Usage

The measured counts must be in the same bins as the prediction of the
[forward model](forward_model.md): same energy range, blinding and binning.
`bin_edges` lists those bins.

```julia
pred = model(params)                 # predicted counts per bin
logL = poisson_likelihood(obs, pred) # obs: measured counts in the same bins
```

- `obs` and `pred` must have the same length, otherwise a `DimensionMismatch` is raised.
- A bin with a predicted count of zero is allowed. A tiny number (`eps()`) is added
  to the prediction inside the logarithm, so the result stays finite.
