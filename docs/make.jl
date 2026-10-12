using Documenter
using Refit

makedocs(;
    sitename = "Refit.jl",
    modules = [Refit],
    pages = [
        "Home" => "index.md",
        "User manual" => [
            "Forward model" => "forward_model.md",
            "Likelihood" => "likelihood.md",
            "Priors" => "priors.md",
        ],
        "API" => "api.md",
    ],
)
