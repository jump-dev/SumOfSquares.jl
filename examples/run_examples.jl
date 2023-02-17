using Test

"""
    _include_sandbox(filename)

Include the `filename` in a temporary module that acts as a sandbox. (Ensuring
no constants or functions leak into other files.)

This function was taken from `JuMP/docs/make.jl`.
"""
function _include_sandbox(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

const _TUTORIAL_DIR = joinpath(@__DIR__, "..", "docs", "src", "tutorials")

@testset "run_examples.jl" begin
    @testset "$dir" for dir in readdir(_TUTORIAL_DIR)
        for filename in readdir(joinpath(_TUTORIAL_DIR, dir))
            if filename[1] == '_'
                continue
            end
            path = joinpath(_TUTORIAL_DIR, dir, filename)
            _include_sandbox(path)
        end
    end
    @testset "Chordal" begin
        include("chordal_sparsity.jl")
        include("chordal_sparsity_with_domain.jl")
    end
end
