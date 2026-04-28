module TestVariableLowRank

using Test
import MultivariateBases as MB
import LowRankOpt as LRO
import Hypatia
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
                            LRO.Factorization(reshape([1.0, 0.0], 2, 1), [1.0]),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(reshape([1.0, 1.0], 2, 1), [1.0]),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(reshape([1.0, 2.0], 2, 1), [1.0]),
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
                            LRO.Factorization(reshape([1.0, 0.0], 2, 1), [2.0]),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(reshape([1.0, 1.0], 2, 1), [2.0]),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization(reshape([1.0, 2.0], 2, 1), [2.0]),
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

function _bridged_model_and_vars()
    @polyvar x
    lag = MB.LagrangeBasis((x,), [[0.0], [1.0], [2.0]])
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    model = MOI.Bridges.Variable.SingleBridgeOptimizer{
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
    }(
        inner,
    )
    set = SumOfSquares.WeightedSOSCone{
        MOI.PositiveSemidefiniteConeTriangle,
    }(
        lag,
        [MB.SubBasis{MB.Monomial}([x^0, x])],
        [MB.algebra_element(1.0 * x^0)],
    )
    p, ci = MOI.add_constrained_variables(model, set)
    return model, inner, p, ci
end

function test_delete()
    model, _, p, _ = _bridged_model_and_vars()
    @test MOI.get(model, MOI.NumberOfVariables()) == 3
    MOI.delete(model, p)
    @test MOI.get(model, MOI.NumberOfVariables()) == 0
    return
end

function test_primal()
    @polyvar x
    lag = MB.LagrangeBasis((x,), [[0.0], [1.0], [2.0]])
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    )
    model = MOI.Bridges.Variable.SingleBridgeOptimizer{
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
    }(
        mock,
    )
    set = SumOfSquares.WeightedSOSCone{
        MOI.PositiveSemidefiniteConeTriangle,
    }(
        lag,
        [MB.SubBasis{MB.Monomial}([x^0, x])],
        [MB.algebra_element(1.0 * x^0)],
    )
    p, ci = MOI.add_constrained_variables(model, set)
    # Mock an optimization result: all inner variables = 1.0
    inner_vars = MOI.get(mock, MOI.ListOfVariableIndices())
    MOI.Utilities.mock_optimize!(
        mock,
        MOI.OPTIMAL,
        (MOI.FEASIBLE_POINT, ones(length(inner_vars))),
    )
    # VariablePrimal through bridge (line 158)
    for vi in p
        @test MOI.get(model, MOI.VariablePrimal(), vi) == 1.0
    end
    # ConstraintPrimal through bridge (line 141)
    cp = MOI.get(model, MOI.ConstraintPrimal(), ci)
    @test length(cp) == 3
    @test all(cp .== 1.0)
    return
end

# Verify that with Hypatia, LagrangeBasis uses SetDotProducts (not PSD)
function test_hypatia_uses_setdotproducts()
    @polyvar x
    optimizer = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        Hypatia.Optimizer{Float64}(),
    )
    MOI.set(optimizer, MOI.Silent(), true)
    model = MOI.Bridges.full_bridge_optimizer(optimizer, Float64)
    MOI.Bridges.add_bridge(
        model,
        SumOfSquares.Bridges.Variable.LowRankBridge{Float64},
    )
    LRO.Bridges.Variable.add_all_bridges(model, Float64)
    lag = MB.LagrangeBasis(
        (x,),
        [[0.0], [1.0], [2.0], [3.0], [4.0]],
    )
    set = SumOfSquares.WeightedSOSCone{
        MOI.PositiveSemidefiniteConeTriangle,
    }(
        lag,
        [MB.SubBasis{MB.Monomial}([x^0, x, x^2])],
        [MB.algebra_element(1.0 * x^0)],
    )
    MOI.add_constrained_variables(model, set)
    cache = optimizer.model_cache
    constraint_types =
        MOI.get(cache, MOI.ListOfConstraintTypesPresent())
    @test any(((F, S),) -> S <: LRO.SetDotProducts, constraint_types)
    @test !any(
        ((F, S),) -> S <: MOI.PositiveSemidefiniteConeTriangle,
        constraint_types,
    )
    return
end

end  # module

TestVariableLowRank.runtests()
