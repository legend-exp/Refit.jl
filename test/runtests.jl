using Aqua
using Refit
include("test_hist_utils.jl")

Aqua.test_all(
    Refit,
    stale_deps = (ignore = [:CairoMakie, :LegendMakie],),
    persistent_tasks = false,
)
