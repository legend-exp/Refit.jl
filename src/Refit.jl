module Refit

using DensityInterface
using IntervalSets
using StatsBase

include("hist_utils.jl")
include("forward_model.jl")
include("likelihood.jl")
include("priors.jl")
include("io.jl")

export rebin_hist, add_histograms, integrate_hist, get_counts, counts_to_histogram
export build_forward_model
export poisson_likelihood, build_likelihood
export get_max_norm, get_priors
export read_hist
end # module Refit
