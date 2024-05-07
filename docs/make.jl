# Inspired from https://github.com/jump-dev/JuMP.jl/blob/master/docs/make.jl

using SumOfSquares
using Documenter, Literate
import DocumenterCitations

const _TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorials")
const _OUTPUT_DIR   = joinpath(@__DIR__, "src", "generated")
const _TUTORIAL_SUBDIR = [
    "Getting started",
    "Polynomial Optimization",
    "Systems and Control",
    "Other Applications",
    "Noncommutative and Hermitian",
    "Sparsity",
    "Symmetry",
    "Extension",
]

function literate_directory(dir)
    output_dir = joinpath(_OUTPUT_DIR, dir)
    for filename in readdir(joinpath(_TUTORIAL_DIR, dir))
        if filename[1] == '_'
            continue
        end
        path = joinpath(_TUTORIAL_DIR, dir, filename)
        Literate.markdown(path, output_dir; documenter = true)
        Literate.notebook(path, output_dir; documenter = true)
        Literate.script(path, output_dir; documenter = true)
    end
end

literate_directory.(_TUTORIAL_SUBDIR)

makedocs(
    sitename = "SumOfSquares",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/citations.css"],
    ),
    pages = [
        "Index" => "index.md",
        "Sum-of-Squares Programming" => "sumofsquares.md",
        "Variables" => "variables.md",
        "Constraints" => "constraints.md",
        "API Reference" => [
            "reference/standard_form.md",
            "reference/constraints.md",
            "reference/certificate.md",
            "reference/internal.md",
        ],
        "Tutorials" => map(
            subdir ->
                subdir => map(
                    file -> joinpath("generated", subdir, file),
                    filter(
                        file -> endswith(file, ".md"),
                        sort(readdir(joinpath(_OUTPUT_DIR, subdir))),
                    ),
                ),
            _TUTORIAL_SUBDIR,
        ),
        "Bibliography" => "bibliography.md",
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [SumOfSquares, PolyJuMP],
    plugins = [
        DocumenterCitations.CitationBibliography(
            joinpath(@__DIR__, "src", "references.bib");
            style = :authoryear,
        ),
    ],
)

deploydocs(
    repo   = "github.com/jump-dev/SumOfSquares.jl.git",
    push_preview = true,
)
