using Aqua
using Refit
using Test

include("test_hist_utils.jl")
include("test_likelihood.jl")

Aqua.test_all(Refit, persistent_tasks = false)
