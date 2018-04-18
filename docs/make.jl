using Documenter, SumOfSquares

makedocs(
    format = :html,
    sitename = "SumOfSquares",
    pages = [
        "Index" => "index.md",
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [SumOfSquares]
)

deploydocs(
    repo   = "github.com/JuliaOpt/SumOfSquares.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
