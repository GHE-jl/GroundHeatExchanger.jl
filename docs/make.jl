using GroundHeatExchanger
using Documenter

# Make `using GroundHeatExchanger` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    GroundHeatExchanger,
    :DocTestSetup,
    :(using GroundHeatExchanger);
    recursive = true,
)

makedocs(;
    modules = [GroundHeatExchanger],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "GroundHeatExchanger.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/GroundHeatExchanger.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Modeling theory" => [
            "Simulation pipeline" => "theory/overview.md",
            "Temporal superposition" => "theory/superposition.md",
            "Fluid temperature" => "theory/temperature.md",
            "g-function compression" => "theory/compression.md",
        ],
        "Utilities" => "utilities.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    # checkdocs only covers symbols defined in GroundHeatExchanger itself; re-exported symbols from
    # BoreholeResistance.jl and GroundResponse.jl are documented in their own package sites.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/GroundHeatExchanger.jl",
    devbranch = "main",
)
