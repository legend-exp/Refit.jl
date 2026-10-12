# Refit.jl

Tools to perform Bayesian fits of binned data with models based on Monte-Carlo
simulations, often generated with [remage](https://remage.readthedocs.io/en/stable/).

`Refit.jl` is a toolkit, it does not provide an end-to-end fitting tool. A fit
compares a measured energy spectrum with a sum of simulated spectra. Each simulated
spectrum is scaled by a free parameter, its activity. The parameters are found by
sampling the posterior, which combines the likelihood of the data with priors.

## User manual

- [Forward model](forward_model.md): map the activities to the predicted counts in each energy bin
- [Likelihood](likelihood.md): compare the predicted counts with the measured ones

## Reference

- [API](api.md): docstrings of all functions
