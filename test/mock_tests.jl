include("Tests/Tests.jl")
include("utilities.jl")

using Test, JuMP

# Needs https://github.com/JuliaOpt/MathOptInterface.jl/pull/669
import Pkg
if Pkg.installed()["MathOptInterface"] > v"0.8.2"
    @testset "Term" begin
        config = MOI.Test.TestConfig()
        optimize!(mock) = MOIU.mock_optimize!(mock, [0.0, 0.0],
            (MOI.VectorAffineFunction{Float64}, MOI.Nonnegatives) => [[1.0]],
            (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0]])
        mock = bridged_mock(optimize!)
        Tests.sos_term_test(mock, config)
        Tests.sdsos_term_test(mock, config)
        Tests.dsos_term_test(mock, config)
    end
end
