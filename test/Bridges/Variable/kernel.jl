module TestVariableKernel

using Test
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

function test_runtests()
    @polyvar x y
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Variable.KernelBridge,
        model -> begin
            MOI.add_constrained_variables(
                model,
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MonomialBasis([x^4, x^3 * y, x^2 * y^2, y^4]),
                    [MonomialBasis([x^2, y^2, x * y])],
                    [1.0 * x^0 * y^0],
                ),
            )
        end,
        model -> begin
            Q, _ = MOI.add_constrained_variables(
                model,
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
        end,
    )
    return
end

end  # module

TestVariableKernel.runtests()
