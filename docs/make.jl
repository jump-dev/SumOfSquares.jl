using SumOfSquares
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

const EXAMPLES = [
    "sos_decomposition.jl",
    "Bound on Global Extremum.jl",
    "Lyapunov Function Search.jl",
    "Sum-of-Squares Matrices.jl",
    "Noncommutative variables.jl",
    "Sums of Hermitian squares.jl",
]

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
end

makedocs(
    sitename = "SumOfSquares",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Index" => "index.md",
        "Sum-of-Squares Programming" => "sumofsquares.md",
        "Variables" => "variables.md",
        "Constraints" => "constraints.md",
        "Examples" => Any[
            "SOS decomposition" => "generated/sos_decomposition.md",
            "Bound on Global Extremum" => "generated/Bound on Global Extremum.md",
            "Lyapunov Function Search" => "generated/Lyapunov Function Search.md",
            "Sum-of-Squares Matrices" => "generated/Sum-of-Squares Matrices.md",
            "Noncommutative variables" => "generated/Noncommutative variables.md",
            "Sums of Hermitian squares" => "generated/Sums of Hermitian squares.md",
        ]
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [SumOfSquares, PolyJuMP]
)

deploydocs(
    repo   = "github.com/jump-dev/SumOfSquares.jl.git",
)
