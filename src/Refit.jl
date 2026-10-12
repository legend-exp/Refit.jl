module Refit

using DensityInterface
using StatsBase

include("hist_utils.jl")
include("forward_model.jl")
include("likelihood.jl")
include("io.jl")

export rebin_hist, add_histograms, integrate_hist, get_counts
export build_forward_model
export poisson_likelihood, build_likelihood
export read_hist
end # module Refit
