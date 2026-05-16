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
