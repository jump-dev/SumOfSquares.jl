module TestVariableKernel

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

# [Parrilo2003, Example 3.5]
function test_runtests()
    @polyvar x y
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Variable.KernelBridge,
        model -> begin
            p, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MB.SubBasis{MB.Monomial}([
                        x^4,
                        x^3 * y,
                        x^2 * y^2,
                        x * y^3,
                        y^4,
                    ]),
                    [MB.SubBasis{MB.Monomial}([x^2, y^2, x * y])],
                    [MB.algebra_element(1.0 * x^0 * y^0)],
                ),
            )
            a = float.(1:length(p))
            MOI.add_constraint(model, MOI.Utilities.vectorize([a' * p]), MOI.Zeros(1))
        end,
        model -> begin
            q, _ = MOI.add_constrained_variables(
                model,
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
            a = float.(1:length(q))
            MOI.add_constraint(model, MOI.Utilities.vectorize([1.0 * q[1] + 2.0 * q[3] + 4.0 * (1.0q[4] + 1.0q[6]) + 6.0 * q[5]]), MOI.Zeros(1))
        end;
        allow_outer_constraint_function_error = true,
    )
    return
end

end  # module

TestVariableKernel.runtests()
