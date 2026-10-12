using Aqua
using Refit
using Test

include("test_hist_utils.jl")
include("test_forward_model.jl")
include("test_likelihood.jl")
include("test_priors.jl")
include("test_workflow.jl")
include("test_io.jl")
Aqua.test_all(Refit, persistent_tasks = false)
