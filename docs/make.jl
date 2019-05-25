using Documenter, SwitchOnSafety

makedocs(
    sitename = "SwitchOnSafety",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/JuliaOpt/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Introduction" => "index.md"
    ]
)

deploydocs(
    repo   = "github.com/blegat/SwitchOnSafety.jl.git"
)
