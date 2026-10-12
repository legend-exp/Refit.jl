using Documenter
using Refit

makedocs(;
    sitename = "Refit.jl",
    modules = [Refit],
    pages = [
        "Home" => "index.md",
        "Forward model" => "forward_model.md",
        "API" => "api.md",
    ],
)
