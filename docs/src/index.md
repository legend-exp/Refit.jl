# Refit.jl

`Refit.jl` contains tools to perform Bayesian (with [BAT.jl](https://bat.github.io/) fits of binned data with models
based on Monte-Carlo simulations, often generated with [remage](https://remage.readthedocs.io/en/stable/).

`Refit.jl` is a toolkit, it does not provide an end-to-end fitting tool.

## User manual
`Refit.jl` can also do a lot more than this, for details see the following pages!

- [Getting started](getting_started.md): install `Refit.jl` and the packages for the fit
- [Reading histograms](reading.md): read the data and the simulated spectra and normalise them
- [Accounting for systematic uncertainties](generalised_hist.md): interpolation of spectra to account for systematic uncertainities in models
- [Forward model](forward_model.md): build models mapping the parameters to the predicted counts in each energy bin
- [Binning](binning.md): tools to handle binning the spectra
- [Likelihood](likelihood.md): functionality to build an optimised likelihood functions
- [Priors](priors.md): including estimation of ranges for the activities from the data
- [Running the fit](fit.md): sample the posterior with BAT.jl and get the best fit
- [Plotting](plotting.md): plot the best fit, make plots with a different style, and plot the posteriors
- [Saving the results](saving.md): save the samples, the best fit spectra and the activities

## Tutorials

- [Fitting a background spectrum](tutorial.md): a complete fit on fake data, from reading the files to sampling the posterior
- [Fitting a peak](tutorial_peak.md): a fit with analytic components only, a Gaussian peak
  on a linear background
- [Simultaneous fits](tutorial_simultaneous.md): two detectors fitted together,
  which removes degeneracies
- [Fitting the FCCD with template morphing](tutorial_fccd.md): the dead layer of a
  detector fitted to a 60 keV γ line, with spectra simulated on a grid

## Reference

- [API](api.md): docstrings of all functions
