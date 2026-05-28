module TestVariableLowRank

using Test
import MultivariateBases as MB
import LowRankOpt as LRO
import Hypatia
import Clarabel
import Dualization
using DynamicPolynomials
using JuMP
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
                            LRO.Factorization([1.0, 0.0], reshape([1.0], ())),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization([1.0, 1.0], reshape([1.0], ())),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization([1.0, 2.0], reshape([1.0], ())),
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
                            LRO.Factorization([1.0, 0.0], reshape([2.0], ())),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization([1.0, 1.0], reshape([2.0], ())),
                        ),
                        LRO.TriangleVectorization(
                            LRO.Factorization([1.0, 2.0], reshape([2.0], ())),
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
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
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
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
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
    lag = MB.LagrangeBasis((x,), [[0.0], [1.0], [2.0], [3.0], [4.0]])
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
        lag,
        [MB.SubBasis{MB.Monomial}([x^0, x, x^2])],
        [MB.algebra_element(1.0 * x^0)],
    )
    MOI.add_constrained_variables(model, set)
    cache = optimizer.model_cache
    constraint_types = MOI.get(cache, MOI.ListOfConstraintTypesPresent())
    @test any(((F, S),) -> S <: LRO.SetDotProducts, constraint_types)
    @test !any(
        ((F, S),) -> S <: MOI.PositiveSemidefiniteConeTriangle,
        constraint_types,
    )
    return
end

# Build the JuMP model used by `test_hypatia_jump_uses_setdotproducts` and
# `test_clarabel_jump_lagrange_via_lro_bridges`: a univariate SOS problem
# `max γ s.t. p - γ ∈ SOS` whose unique optimum is `γ = -6`, with a
# `BoxSampling` `zero_basis` so the basis carried by `WeightedSOSCone` is a
# `LagrangeBasis` and `Variable.LowRankBridge` becomes the relevant entry
# point. Using a shared helper guarantees the two solvers receive the exact
# same model.
function _lagrange_jump_model(
    optimizer;
    optimize = true,
    pre_constraint = nothing,
)
    @polyvar x
    model = JuMP.Model(optimizer)
    JuMP.set_silent(model)
    if pre_constraint !== nothing
        pre_constraint(model)
    end
    @variable(model, γ)
    @objective(model, Max, γ)
    p = x^4 - 4x^3 - 2x^2 + 12x + 3
    @constraint(
        model,
        p - γ in SOSCone(),
        zero_basis = MB.BoxSampling([-1.0], [1.0]),
    )
    if optimize
        JuMP.optimize!(model)
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.value(γ) ≈ -6 rtol = 1e-4
    else
        # No sub-solver loaded; just push the constraints down the bridge layer
        # so we can inspect what the inner optimizer received.
        MOI.Utilities.attach_optimizer(JuMP.backend(model))
    end
    inner_cache = JuMP.backend(model).optimizer.model.model_cache
    return MOI.get(inner_cache, MOI.ListOfConstraintTypesPresent())
end

# Driven through the JuMP `@constraint` macro and a `BoxSampling` `zero_basis`,
# the SOS constraint should reach Hypatia as a `VOV-in-LRO.SetDotProducts`
# rather than going through `MOI.PositiveSemidefiniteConeTriangle` (which would
# require Hypatia to bridge it back to its native scaled cone).
#
# The bridge chain that lands the constraint at Hypatia is:
#   1. `PolyJuMP`'s `bridges(::Type{<:WeightedSOSCone})` auto-registers
#      `SumOfSquares.Bridges.Variable.LowRankBridge`, which produces
#      `LRO.SetDotProducts{WITHOUT_SET, PSDConeTriangle,
#       TriangleVectorization{T, Factorization{T, Vector{T}, Array{T,0}}}}` —
#      i.e. directly emits per-Lagrange-point rank-1 factors `U[j, :]` with a
#      0-dim scaling `weights[j]` (rather than baking the weight into the
#      vector as `sqrt(weights[j]) * column`, which would force a `sqrt`
#      that's only valid for non-negative weights).
#   2. `PolyJuMP.bridges` for that rank-1 `SetDotProducts` variant adds
#      `LRO.Bridges.Variable.ToPositiveBridge`, which collapses the 0-dim
#      scaling into `LRO.One{T}` (= `FillArrays.Ones{T,0,Tuple{}}`).
#   3. Hypatia's MOI wrapper natively consumes that final rank-1 variant.
function test_hypatia_jump_uses_setdotproducts()
    constraint_types = _lagrange_jump_model(Hypatia.Optimizer)
    @test !any(
        ((F, S),) -> S <: MOI.PositiveSemidefiniteConeTriangle,
        constraint_types,
    )
    # The constraint reaches Hypatia as the rank-1 `LRO.SetDotProducts` variant.
    # Hypatia's MOI wrapper consumes this set as a `VectorAffineFunction`
    # constraint (the `VectorOfVariables` form is mapped to `VAF` by
    # `MOI.Bridges.Constraint.VectorFunctionizeBridge`). The point of this
    # assertion is that the rank-1 low-rank structure survives the bridge
    # chain — it is *not* collapsed back to a dense `PSDConeTriangle`.
    expected_S = LRO.SetDotProducts{
        LRO.WITHOUT_SET,
        MOI.PositiveSemidefiniteConeTriangle,
        LRO.TriangleVectorization{
            Float64,
            LRO.Factorization{Float64,Vector{Float64},LRO.One{Float64}},
        },
    }
    @test (MOI.VectorAffineFunction{Float64}, expected_S) in constraint_types
    return
end

# Counterpart to `test_hypatia_jump_uses_setdotproducts`: for a classical SDP
# solver (Clarabel) without a low-rank interface, `KernelBridge` is not
# applicable for the `LagrangeBasis` case (its `supports_constrained_variable`
# returns `false`). The bridge graph instead routes `LowRankBridge`'s LRO
# output back to classical PSD through the `AppendSetBridge` →
# `DotProductsBridge` fallback registered in `src/variables.jl`. Clarabel
# itself supports the scaled PSD cone, so MOI's own `SetDotScalingBridge`
# finishes the chain.
function test_clarabel_jump_lagrange_via_lro_bridges()
    constraint_types = _lagrange_jump_model(Clarabel.Optimizer)
    # The LRO chain (`AppendSetBridge` + `DotProductsBridge`) lowers
    # `LRO.SetDotProducts` back to a classical PSD cone, which MOI's
    # `SetDotScalingBridge` then maps to Clarabel's natively-supported
    # `MOI.Scaled{PSDConeTriangle}`. No `LRO.SetDotProducts` reaches Clarabel.
    @test !any(((F, S),) -> S <: LRO.SetDotProducts, constraint_types)
    @test (
        MOI.VectorAffineFunction{Float64},
        MOI.Scaled{MOI.PositiveSemidefiniteConeTriangle},
    ) in constraint_types
    return
end

# Same `_lagrange_jump_model` setup but with `LRO.Optimizer` as the inner
# optimizer. `LRO.Optimizer` supports `VAF`-in-`LinearCombinationInSet` (its
# low-rank constraint form) *and* `VAF`-in-`PSDConeTriangle` (its classical
# form), so without intervention the bridge layer prefers the cheaper one-step
# `ImageBridge` path to classical PSD and the rank-1 structure carried by
# `Variable.LowRankBridge` never reaches the solver.
#
# Workaround: wrap the inner optimizer in `Dualization.dual_optimizer` and
# remove `Constraint.ImageBridge` from the outer bridge layer. With the
# constraint-side `ImageBridge` shortcut gone, the cost graph routes the SOS
# constraint through `Variable.LowRankBridge` → `LRO.SetDotProducts`. The
# `dual_optimizer` wrapper is what makes the LRO solver actually see this
# low-rank-structured cone: `LRO.Optimizer` natively consumes
# `LinearCombinationInSet` (constraint-side LRO), so the dualization layer is
# needed to flip the `SetDotProducts` (variable-side LRO) into its constraint
# dual. The `remove_bridge` shim is documented in
# `docs/src/tutorials/Noncommutative and Hermitian/chsh.jl` and tracked by
# <https://github.com/jump-dev/MathOptInterface.jl/pull/3001>; once that lands
# we can drop the `remove_bridge` call. The `dual_optimizer` wrapper remains
# required regardless, since it bridges the variable/constraint divide that
# the LRO bridge graph does not itself cross.
function test_lro_optimizer_preserves_low_rank_structure()
    constraint_types = _lagrange_jump_model(
        Dualization.dual_optimizer(LRO.Optimizer{Float64});
        optimize = false,
        pre_constraint = model -> begin
            backend = JuMP.backend(model)
            SumOfSquares.Bridges.add_all_bridges(backend.optimizer, Float64)
            MOI.Bridges.remove_bridge(
                backend.optimizer,
                SumOfSquares.Bridges.Constraint.ImageBridge{Float64},
            )
        end,
    )
    @test any(((F, S),) -> S <: LRO.SetDotProducts, constraint_types)
    @test !any(
        ((F, S),) -> S <: MOI.PositiveSemidefiniteConeTriangle,
        constraint_types,
    )
    return
end

end  # module

TestVariableLowRank.runtests()
