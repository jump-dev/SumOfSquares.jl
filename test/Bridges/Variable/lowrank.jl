module TestVariableLowRank

using Test
import MultivariateBases as MB
import LowRankOpt as LRO
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

# 1 variable, 3 Lagrange points, gram basis [1, x], weight 1
# U = transformation_to([1, x], LagrangeBasis at [0, 1, 2]) = [1 0; 1 1; 1 2]
# weights at points = [1, 1, 1]
function test_runtests()
    @polyvar x
    lag = MB.LagrangeBasis((x,), [[0.0], [1.0], [2.0]])
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Variable.LowRankBridge,
        model -> begin
            p, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    lag,
                    [MB.SubBasis{MB.Monomial}([x^0, x])],
                    [MB.algebra_element(1.0 * x^0)],
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
                LRO.SetDotProducts{LRO.WITHOUT_SET}(
                    SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
                    [
                        LRO.TriangleVectorization(
                            LRO.Factorization(
                                reshape([1.0, 0.0], 2, 1),
                                [1.0],
                            ),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(
                                reshape([1.0, 1.0], 2, 1),
                                [1.0],
                            ),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(
                                reshape([1.0, 2.0], 2, 1),
                                [1.0],
                            ),
                        ),
                    ],
                ),
            )
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([1.0 * q[1] + 2.0 * q[2] + 3.0 * q[3]]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# Same as above but with weight 2
function test_runtests_weighted()
    @polyvar x
    lag = MB.LagrangeBasis((x,), [[0.0], [1.0], [2.0]])
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Variable.LowRankBridge,
        model -> begin
            p, _ = MOI.add_constrained_variables(
                model,
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    lag,
                    [MB.SubBasis{MB.Monomial}([x^0, x])],
                    [MB.algebra_element(2.0 * x^0)],
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
                LRO.SetDotProducts{LRO.WITHOUT_SET}(
                    SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
                    [
                        LRO.TriangleVectorization(
                            LRO.Factorization(
                                reshape([1.0, 0.0], 2, 1),
                                [2.0],
                            ),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(
                                reshape([1.0, 1.0], 2, 1),
                                [2.0],
                            ),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(
                                reshape([1.0, 2.0], 2, 1),
                                [2.0],
                            ),
                        ),
                    ],
                ),
            )
            MOI.add_constraint(
                model,
                MOI.Utilities.vectorize([1.0 * q[1] + 2.0 * q[2] + 3.0 * q[3]]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

end  # module

TestVariableLowRank.runtests()
