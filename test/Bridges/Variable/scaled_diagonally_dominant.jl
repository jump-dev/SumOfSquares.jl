module TestVariableScaledDiagonallyDominant

using Test
import MultivariateBases as MB
using DynamicPolynomials
using SumOfSquares

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_error_dim_1()
    model = MOI.Utilities.Model{Float64}()
    bridged = MOI.Bridges.Variable.SingleBridgeOptimizer{
        SumOfSquares.Bridges.Variable.ScaledDiagonallyDominantBridge{Float64},
    }(
        model,
    )
    err = ErrorException(
        "The bridges does not work with 1, `matrix_cone` should have returned `Nonnegatives` instead.",
    )
    @test_throws err MOI.add_constrained_variables(
        bridged,
        SumOfSquares.ScaledDiagonallyDominantConeTriangle(1),
    )
end

function test_runtests()
    @polyvar x y
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Variable.ScaledDiagonallyDominantBridge,
        model -> begin
            p, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.ScaledDiagonallyDominantConeTriangle(3),
            )
        end,
        model -> begin
            q1, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
            )
            q2, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
            )
            q3, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
            )
        end,
    )
    return
end

end  # module

TestVariableScaledDiagonallyDominant.runtests()
