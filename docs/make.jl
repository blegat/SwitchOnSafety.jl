using SwitchOnSafety
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

const EXAMPLES = [
    "AJPR14e54.jl",
]

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
end

makedocs(
    sitename = "SwitchOnSafety",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Index" => "index.md",
        "Invariant Sets" => "invariant.md",
        "Joint Spectral Radius" => "jsr.md",
        "Examples" => Any[
            "AJPR14e54" => "generated/AJPR14e54.md",
        ]
    ]
)

deploydocs(
    repo   = "github.com/blegat/SwitchOnSafety.jl.git",
    push_preview = true,
)
