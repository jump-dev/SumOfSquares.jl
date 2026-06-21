module TestConstraintImage

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

# [Parrilo2003, Example 6.1]
function test_runtests()
    T = Float64
    @polyvar x y
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.ImageBridge,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(
                    MOI.VectorAffineTerm{T}[],
                    T[5, -1, 2, 2],
                ),
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MB.SubBasis{MB.Monomial}([y^4, x^2 * y^2, x^3 * y, x^4]),
                    [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
                    [MB.algebra_element(one(T) * x^0 * y^0)],
                ),
            )
        end,
        model -> begin
            λ = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.Utilities.operate(
                    vcat,
                    T,
                    T(5),
                    T(0),
                    T(2) * λ - T(1),
                    T(-1) * λ,
                    T(1),
                    T(2),
                ),
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
        end,
    )
    return
end

# Same as `test_runtests` but with a non-unit constant weight. This is the
# kind of `WeightedSOSCone` that arises when the `domain` keyword is used in
# `@constraint` and the resulting Putinar-style decomposition leaves a scaling
# in front of one of the SOS multipliers. Mirrors
# `test_runtests_weighted` in `test/Bridges/Variable/kernel.jl`.
function test_runtests_weighted()
    T = Float64
    @polyvar x y
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.ImageBridge,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(
                    MOI.VectorAffineTerm{T}[],
                    T[5, -1, 2, 2],
                ),
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MB.SubBasis{MB.Monomial}([y^4, x^2 * y^2, x^3 * y, x^4]),
                    [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
                    [MB.algebra_element(T(2) * x^0 * y^0)],
                ),
            )
        end,
        model -> begin
            λ = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.Utilities.operate(
                    vcat,
                    T,
                    T(5 // 2),
                    T(0),
                    T(2) * λ - T(1 // 2),
                    T(-1) * λ,
                    T(1 // 2),
                    T(1),
                ),
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# Same as `test_runtests` but with a non-constant weight, which is the form
# that arises when the `domain` keyword introduces a polynomial multiplier
# `g(x)` (e.g. coming from `@set g >= 0`) into the SOS certificate.
# Mirrors `test_runtests_polynomial_weight` in
# `test/Bridges/Variable/kernel.jl`.
#
# weight (1 + x), gram basis [1, x]
# (1 + x) * (Q11 + 2*Q12*x + Q22*x^2)
#   = Q11 + (2*Q12 + Q11)*x + (Q22 + 2*Q12)*x^2 + Q22*x^3
# so the gram entries are determined (modulo a single zero constraint) by:
#   Q11   = scalars[1]
#   2*Q12 = scalars[2] - scalars[1]
#   Q22   = scalars[3] - scalars[2] + scalars[1]
# with the consistency requirement
#   scalars[1] - scalars[2] + scalars[3] - scalars[4] = 0.
function test_runtests_polynomial_weight()
    T = Float64
    @polyvar x
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.ImageBridge,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(
                    MOI.VectorAffineTerm{T}[],
                    T[1, 3, 4, 2],
                ),
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MB.SubBasis{MB.Monomial}([x^0, x, x^2, x^3]),
                    [MB.SubBasis{MB.Monomial}([x^0, x])],
                    [MB.algebra_element(T(1) * x^0 + T(1) * x)],
                ),
            )
        end,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(
                    MOI.VectorAffineTerm{T}[],
                    T[1, 1, 2],
                ),
                SumOfSquares.PositiveSemidefinite2x2ConeTriangle(),
            )
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[0]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# Same polynomial weight as `test_runtests_polynomial_weight` but with a
# 3-element gram basis `[1, x, x^2]`. Here the greedy contribution order
# inside an entry is sub-optimal: for entry `(1, 3)` the algorithm sees
# `x^2` first (already anchored by entry `(2, 2)`) and so introduces a slack
# `λ`, even though processing `(1, 3)`'s second contribution (`x^3`,
# currently unanchored) first would have allowed it to anchor a fresh
# monomial and use the `entry_anchored = true` subtract branch for `x^2`.
# Bipartite max-cardinality matching between gram entries and basis
# monomials would yield 6 anchors and zero slacks here; the current greedy
# yields 5 anchors and 1 slack, with `x^5` becoming a zero constraint.
function test_runtests_polynomial_weight_avoidable_slack()
    T = Float64
    @polyvar x
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.ImageBridge,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(
                    MOI.VectorAffineTerm{T}[],
                    T[1, 3, 5, 5, 3, 1],
                ),
                SumOfSquares.WeightedSOSCone{
                    MOI.PositiveSemidefiniteConeTriangle,
                }(
                    MB.SubBasis{MB.Monomial}([x^0, x, x^2, x^3, x^4, x^5]),
                    [MB.SubBasis{MB.Monomial}([x^0, x, x^2])],
                    [MB.algebra_element(T(1) * x^0 + T(1) * x)],
                ),
            )
        end,
        model -> begin
            λ = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.Utilities.operate(
                    vcat,
                    T,
                    T(1),
                    T(1),
                    T(3) + T(2) * λ,
                    T(-1) * λ,
                    T(1),
                    T(1),
                ),
                MOI.PositiveSemidefiniteConeTriangle(3),
            )
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[0]),
                MOI.Zeros(1),
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# Regression test for an oversight in
# `MOI.Bridges.added_constraint_types(ImageBridge)`: the bridge can produce a
# constraint in any of the four sets returned by `SOS.matrix_cone(M, k)` for
# `k = 0, 1, 2, 3` depending on the gram-basis length, but only
# `PositiveSemidefiniteConeTriangle` (i.e. `k = 3`) was advertised. The other
# three variants were missing, leaving `LazyBridgeOptimizer` with an
# inconsistent view of the bridge graph (in particular, its bridging cost was
# under-estimated relative to `Variable.KernelBridge`, which advertises all
# four).
function test_added_constraint_types_covers_all_matrix_cones()
    T = Float64
    F = MOI.VectorAffineFunction{T}
    G = MOI.VectorAffineFunction{T}
    M = MOI.PositiveSemidefiniteConeTriangle
    BT = SumOfSquares.Bridges.Constraint.ImageBridge{T,F,G,M}
    added = MOI.Bridges.added_constraint_types(BT)
    for k in 0:3
        @test (F, typeof(SumOfSquares.matrix_cone(M, k))) in added
    end
    @test (G, MOI.Zeros) in added
    return
end

end  # module

TestConstraintImage.runtests()
