# Refit.jl

[![Codecov](https://codecov.io/gh/legend-exp/Refit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/legend-exp/Refit.jl)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)


Tools to perform Bayesian fits of binned data with models based on Monte-Carlo simulations, often generated with [remage](https://remage.readthedocs.io/en/stable/).

This includes:
- Reading and manipulation of Monte-Carlo files,
- Automatic handeling of priors and initial parameter ranges,
- Dedicated methods to incorportate systematic uncertainties,
- Plotting and saving outputs.


> Warning: `Refit.jl` is still under heavy development, use at your risk!
