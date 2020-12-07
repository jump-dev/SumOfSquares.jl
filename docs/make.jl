using SumOfSquares
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

include(joinpath(EXAMPLES_DIR, "all_examples.jl"))
deleteat!(EXAMPLES, findfirst(isequal("chordal_sparsity.jl"), EXAMPLES))
deleteat!(EXAMPLES, findfirst(isequal("chordal_sparsity_with_domain.jl"), EXAMPLES))

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
        "Examples" => map(EXAMPLES) do jl_file
            # Need `string` as Documenter fails if `name` is a `SubString{String}`.
            name = string(split(jl_file, ".")[1])
            return name => "generated/$name.md"
        end
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [SumOfSquares, PolyJuMP]
)

deploydocs(
    repo   = "github.com/jump-dev/SumOfSquares.jl.git",
    push_preview = true,
)
