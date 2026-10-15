# no display is needed to draw the plots of the tutorial
ENV["GKSwstype"] = "100"

using Documenter
using Refit

makedocs(;
    sitename = "Refit.jl",
    format = Documenter.HTML(; canonical = "https://legend-exp.github.io/Refit.jl"),
    modules = [Refit],
    pages = [
        "Home" => "index.md",
        "User manual" => [
            "Reading spectra" => "reading.md",
            "Forward model" => "forward_model.md",
            "Binning" => "binning.md",
            "Likelihood" => "likelihood.md",
            "Priors" => "priors.md",
            "Plotting" => "plotting.md",
            "Saving the results" => "saving.md",
        ],
        "Tutorials" => [
            "Fitting a background spectrum" => "tutorial.md",
            "Fitting a peak" => "tutorial_peak.md",
            "Simultaneous fits" => "tutorial_simultaneous.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/legend-exp/Refit.jl.git",
    devbranch = "main",
    push_preview = true,
)
