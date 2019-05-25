using Documenter, SwitchOnSafety

makedocs(
    sitename = "SwitchOnSafety",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Index" => "index.md",
        "Invariant Sets" => "invariant.md",
        "Joint Spectral Radius" => "jsr.md"
    ]
)

deploydocs(
    repo   = "github.com/blegat/SwitchOnSafety.jl.git"
)
