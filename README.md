# Refit.jl

[![Codecov](https://codecov.io/gh/legend-exp/Refit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/legend-exp/Refit.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://legend-exp.github.io/Refit.jl/dev)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)


Tools to perform Bayesian fits of binned data with models based on Monte-Carlo simulations, often generated with [remage](https://remage.readthedocs.io/en/stable/).

This includes:
- Reading and manipulation of Monte-Carlo files,
- Automatic handeling of priors and initial parameter ranges,
- Dedicated methods to incorportate systematic uncertainties,
- Plotting and saving outputs.


> [!WARNING]
> `Refit.jl` is still under heavy development, use at your risk!


## Basic usage
`Refit.jl` is a toolkit, it does not provide an end-to-end fitting tool, but instead a set of functionality.

A minimal fit can be performed with:

```julia

using Refit

# read the data as a Histogram
data = read_hist("hist/spectrum",data_file)

# dictionary of MC pdfs to fit, the keys are Symbols and name the parameters
mc_dict = Dict(
    :Ar39_liquid_argon => read_pdf("hist/spectrum","$mc_files/Ar39_liquid_argon.lh5"),
    :two_vbb_hpge_bulk => read_pdf("hist/spectrum","$mc_files/2vbb_hpge_bulk.lh5"),
    # ...
)

# define priors (by hand, from input files or from the data), with the same keys
priors = Dict()
priors[:Ar39_liquid_argon] = 0.. 2000
priors[:two_vbb_hpge_bulk] = 0.. 2
# ...

# or guess ranges for the fit from the data
# priors = get_priors(mc_dict, data; e_low = 500, e_high = 3200)

# convert to a distribution
prior_dist = build_prior_dist(priors)

# get a model mapping parameters to data, and the edges of the bins it predicts
model, bin_edges = build_forward_model(mc_dict; bin = 20, e_low = 500, e_high = 3200)

# build a likelihood / posterior
counts = get_counts(data, bin_edges)
likelihood = build_likelihood(counts, model)
posterior = PosteriorMeasure(likelihood, prior_dist)

# sample
samples =  bat_sample(....)

# plots saving etc
```

The forward model is described in the [user manual](https://legend-exp.github.io/Refit.jl/dev/).

More advanced functionality includes:
- Extraction of reasonable uniform priors from the data
- Variation of simulated Histograms based on systematic uncertainties
- Simulataneous fits of multiple data spectra
