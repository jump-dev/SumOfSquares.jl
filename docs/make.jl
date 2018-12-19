using Documenter, SumOfSquares

makedocs(
    sitename = "SumOfSquares",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    html_prettyurls = get(ENV, "CI", nothing) == "true",
    pages = [
        "Index" => "index.md",
        "Sum-of-Squares Programming" => "sumofsquares.md",
        "Variables" => "variables.md",
        "Constraints" => "constraints.md",
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [SumOfSquares]
)

deploydocs(
    repo   = "github.com/JuliaOpt/SumOfSquares.jl.git",
)
