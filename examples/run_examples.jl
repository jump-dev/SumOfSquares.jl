using Test

const _TUTORIAL_DIR = joinpath(@__DIR__, "..", "docs", "src", "tutorials")

@testset "run_examples.jl" begin
    @testset "$dir" for dir in readdir(_TUTORIAL_DIR)
        for filename in readdir(joinpath(_TUTORIAL_DIR, dir))
            if filename[1] == '_'
                continue
            end
            path = joinpath(_TUTORIAL_DIR, dir, filename)
            include(path)
        end
    end
    @testset "Chordal" begin
        include("chordal_sparsity.jl")
        include("chordal_sparsity_with_domain.jl")
    end
end
