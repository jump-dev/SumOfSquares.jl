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
@testset "Bivariate quadratic" begin
    config = MOI.Test.TestConfig()
    optimize!(mock) = MOIU.mock_optimize!(mock, [2.0, 1.0, 1.0, 1.0],
        (MOI.VectorAffineFunction{Float64}, MOI.RotatedSecondOrderCone) => [[1.0, 1.0, -√2]],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[1.0, -1.0, 1.0]])
    mock = bridged_mock(optimize!)
    Tests.sos_bivariate_quadratic_test(mock, config)
    Tests.sdsos_bivariate_quadratic_test(mock, config)
    optimize!(mock) = MOIU.mock_optimize!(mock, [2.0, 1.0, 1.0, 1.0, 1.0],
        (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}) => [0.0, 2.0, 2.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[2.0, -1.0, 0.0]])
    Tests.dsos_bivariate_quadratic_test(mock, config)
end
