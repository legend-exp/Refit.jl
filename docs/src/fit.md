# Running the fit

`Refit.jl` builds the [likelihood](likelihood.md) and the [prior](priors.md).
[BAT.jl](https://bat.github.io/BAT.jl/stable/) samples the posterior, which combines
the two. The main commands are below, with the model, data and priors of the other
examples in this manual:

```@setup fit
include("example_data.jl")
using IntervalSets

model, bin_edges =
    build_forward_model(mc_dict; e_low = 500, e_high = 3200, norm = 1e-3 * livetime)
likelihood = build_likelihood(get_counts(data, bin_edges), model)
ranges = get_priors(mc_dict, data; e_low = 500, e_high = 3200, norm = 1e-3 * livetime)
prior = build_prior_dist(ranges)
```

## Sampling the posterior

`bat_sample` draws samples from the posterior with a Markov chain Monte Carlo (MCMC):
`nchains` chains of `nsteps` steps each.

```@example fit
using BAT

posterior = PosteriorMeasure(likelihood, prior)
samples = bat_sample(
    posterior,
    TransformedMCMC(proposal = RandomWalk(), nsteps = 10^4, nchains = 4),
).result
nothing # hide
```

A real fit needs more steps, see the [tutorials](tutorial.md) for typical settings. The
[BAT.jl documentation](https://bat.github.io/BAT.jl/stable/) describes the sampling
algorithms, the checks of convergence and the output.

## Best fit

The best fit is the global mode of the samples, the sample with the largest posterior
density:

```@example fit
best = mode(samples)
```

The mode can also be found without sampling, by maximising the posterior with an
optimiser. `bat_findmode` needs the package `Optim` for this:

```julia
using Optim

best = bat_findmode(posterior, OptimAlg(optalg = Optim.NelderMead())).result
```

The uncertainty of each parameter, the smallest interval that contains 68% of its
posterior, comes from the samples (see [Saving the results](saving.md)).
