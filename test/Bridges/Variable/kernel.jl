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
                        y^4,
                        x * y^3,
                        x^2 * y^2,
                        x^3 * y,
                        x^4,
                    ]),
                    [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
                    [MB.algebra_element(1.0 * x^0 * y^0)],
                ),
            )
            a = float.(1:length(p))
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([a' * p]),
                MOI.Zeros(1),
            )
        end,
        model -> begin
            q, _ = MOI.add_constrained_variables(
                model,
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
            a = float.(1:length(q))
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([
                    1.0 * q[1] +
                    4.0 * q[2] +
                    3.0 * (1.0q[3] + 2.0q[4]) +
                    8.0 * q[5] +
                    5.0 * q[6],
                ]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

function test_runtests_weighted()
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
                        y^4,
                        x * y^3,
                        x^2 * y^2,
                        x^3 * y,
                        x^4,
                    ]),
                    [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
                    [MB.algebra_element(2.0 * x^0 * y^0)],
                ),
            )
            a = float.(1:length(p))
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([a' * p]),
                MOI.Zeros(1),
            )
        end,
        model -> begin
            q, _ = MOI.add_constrained_variables(
                model,
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
            a = float.(1:length(q))
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([
                    2.0 * q[1] +
                    8.0 * q[2] +
                    6.0 * (1.0q[3] + 2.0q[4]) +
                    16.0 * q[5] +
                    10.0 * q[6],
                ]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# weight (1 + x), gram basis [1, x]
# (1 + x) * (Q11 + 2*Q12*x + Q22*x^2)
#   = Q11 + (2*Q12 + Q11)*x + (Q22 + 2*Q12)*x^2 + Q22*x^3
function test_runtests_polynomial_weight()
    @polyvar x
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Variable.KernelBridge,
        model -> begin
            p, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MB.SubBasis{MB.Monomial}([x^0, x, x^2, x^3]),
                    [MB.SubBasis{MB.Monomial}([x^0, x])],
                    [MB.algebra_element(1.0 * x^0 + 1.0 * x)],
                ),
            )
            a = float.(1:length(p))
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([a' * p]),
                MOI.Zeros(1),
            )
        end,
        model -> begin
            q, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
            )
            # p[1] = 1.0*q[1], p[2] = 1.0*q[1] + 2.0*q[2],
            # p[3] = 2.0*q[2] + 1.0*q[3], p[4] = 1.0*q[3]
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([
                    1.0 * q[1] +
                    2.0 * (1.0q[1] + 2.0q[2]) +
                    3.0 * (2.0q[2] + 1.0q[3]) +
                    4.0 * q[3],
                ]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

end  # module

TestVariableKernel.runtests()
