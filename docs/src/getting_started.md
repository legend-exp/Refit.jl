# Getting started

## Installation

`Refit.jl` needs Julia 1.12 or newer. It is not in the Julia General registry, so it
is installed from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/legend-exp/Refit.jl")
```

We suggest to install it in a project environment (a folder with its own list of packages and
versions) rather than in the global one. This keeps the versions of your analysis
fixed and separate from other work:

```julia
Pkg.activate(".")
Pkg.add(url = "https://github.com/legend-exp/Refit.jl")
```

## Packages for the fit

`Refit.jl` builds the likelihood and the priors, but it does not sample the posterior
or find the best fit. These steps use other packages, which must be installed
separately:

```julia
Pkg.add(name = "BAT", version = "5")
Pkg.add(["Optim", "ForwardDiff"])
```

- [BAT.jl](https://bat.github.io/) samples the posterior. The examples in this manual
  need BAT 5.
- [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) finds the best fit
  (the mode of the posterior).
- [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl/stable/) computes the
  derivatives of the likelihood with respect to the parameters, for optimisers
  that use the gradient, such as LBFGS.

They are not dependencies of `Refit.jl` because the likelihood follows the
`DensityInterface` standard, so any sampler or optimiser can use it. You choose the
packages and their versions.

## First steps

The [tutorial](tutorial.md) runs a complete fit, from reading the files to plotting
the results. The pages of this user manual explain each step in detail.
