# Running the fit

`Refit.jl` builds the [likelihood](likelihood.md) and the [prior](priors.md).
[BAT.jl](https://bat.github.io/BAT.jl/stable/) samples the posterior, which combines
the two. The examples in this manual use BAT.jl version 5. The main commands are below,
with the model, data and priors of the other examples in this manual:

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
optimiser. `bat_findmode` needs the package `Optim` for this. The models and likelihoods
of `Refit.jl` work with automatic differentiation by `ForwardDiff`, so an optimiser that
uses the gradient, such as LBFGS, can be used:

```@example fit
using Optim, ForwardDiff

context = BATContext(ad = ForwardDiff)
best_lbfgs = bat_findmode(posterior, OptimAlg(optalg = Optim.LBFGS()), context).result
```

!!! note "Maximum likelihood"
    `bat_bgml` finds the maximum of the likelihood instead of the posterior. The prior
    only chooses the parameter space of the search, so the result does not depend on it
    as long as it is not zero where the likelihood is large:

    ```@example fit
    best_ml = bat_bgml(likelihood, prior, OptimAlg(optalg = Optim.LBFGS()), context).result
    ```

The uncertainty of each parameter, the smallest interval that contains 68% of its
posterior, comes from the samples (see [Saving the results](saving.md)).

## Quantities derived from the samples

Any quantity that follows from the predicted counts has a posterior: its values for the
samples. `get_predicted_counts(samples, model, bin_edges, ranges)` gives the counts
predicted in energy ranges, in keV, for `n` samples (default 10000) drawn with a
probability proportional to their weight. The result has the total in `:total`, followed
by each source. The ends of the ranges must be edges of the bins of the model. For
example, the counts in the ``^{40}``K line at 1461 keV:

```@example fit
using Statistics

counts = get_predicted_counts(samples, model, bin_edges, 1450 .. 1470)
quantile(counts[:total], [0.16, 0.5, 0.84])  # median and central 68% interval
```

`group_histograms(counts, groups)` sums the sources into groups, sample by sample, and
[`plot_posteriors!`](@ref) draws the posteriors (see [Plotting](plotting.md)).

`background_index(samples, model, bin_edges; livetime, mass)` gives the posterior of the
background index at the Q-value of the 0νββ decay of ``^{76}``Ge, in counts/(keV kg yr):
the counts in a window around 2039 keV, divided by its width and by the exposure, the
product of `livetime` and `mass`, given with units. The default window goes from 1930 to
2190 keV without the γ lines at 2104 and 2119 keV, so the bins of the model must have
edges at 1930, 2099, 2109, 2114, 2124 and 2190 keV, for example with
`binning = [[1930, 2099], [2109, 2114], [2124, 2190]]`.

Samples drawn from the prior instead of the posterior give these quantities before a fit,
for example the background expected from the measured activities of the components (see
[Priors from measured activities](priors.md#Priors-from-measured-activities)):

```julia
samples = (v = [rand(prior) for _ in 1:10^5], weight = ones(10^5))
```
