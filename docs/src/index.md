# Refit.jl

Tools to perform Bayesian (with [BAT.jl](https://bat.github.io/) fits of binned data with models
based on Monte-Carlo simulations, often generated with [remage](https://remage.readthedocs.io/en/stable/).

`Refit.jl` is a toolkit, it does not provide an end-to-end fitting tool. A fit
compares a measured energy spectrum with a sum of simulated spectra. Each simulated
spectrum is scaled by the activity. The parameters are found by
sampling the posterior, which combines the likelihood of the data with priors.

## User manual
`Refit.jl` can also do a lot more than this, for details see the following pages!

- [Reading histograms](reading.md): read the data and the simulated spectra and normalise them
- [Forward model](forward_model.md): map the activities to the predicted counts in each energy bin
- [Likelihood](likelihood.md): compare the predicted counts with the measured ones
- [Priors](priors.md): estimate ranges for the activities from the data
- [Plotting](plotting.md): plot the best fit, make plots with a different style, and plot the posteriors
- [Saving the results](saving.md): save the samples, the best fit spectra and the activities

## Tutorials

- [Fitting a background spectrum](tutorial.md): a complete fit on made-up data, from reading the files to sampling the posterior
- [Fitting a peak](tutorial_peak.md): a fit with analytic components only, a Gaussian peak
  on a linear background
- [Simultaneous fits](tutorial_simultaneous.md): two detectors fitted together,
  which separates two sources that one detector alone cannot tell apart

## Reference

- [API](api.md): docstrings of all functions
