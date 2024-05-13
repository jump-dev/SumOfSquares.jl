using Test

@testset "$(dir)" for dir in ["Variable", "Constraint"]
    @testset "$(file)" for file in readdir(joinpath(@__DIR__, dir))
        if !endswith(file, ".jl")
            continue
        elseif dir == "." && file == "Bridges.jl"
            continue
        end
        include(joinpath(@__DIR__, dir, file))
    end
end
