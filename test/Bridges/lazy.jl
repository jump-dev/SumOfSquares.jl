module TestLazy

using Test
import Clarabel
using DynamicPolynomials
import MultivariateBases as MB
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

# Verify that with a full bridge optimizer wrapping Clarabel, the bridges added by
# `SumOfSquares.Bridges.add_all_bridges` route a `WeightedSOSCone` constraint
# through `ImageBridge` (Clarabel natively supports `MOI.PositiveSemidefiniteConeTriangle`
# as a constraint, so the image-form bridge is the cheapest path).
function test_clarabel_uses_image_bridge()
    T = Float64
    F = MOI.VectorAffineFunction{T}
    @polyvar x y
    optimizer = MOI.instantiate(
        Clarabel.Optimizer;
        with_cache_type = T,
        with_bridge_type = T,
    )
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
        MB.SubBasis{MB.Monomial}([y^4, x^2 * y^2, x^3 * y, x^4]),
        [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
        [MB.algebra_element(one(T) * x^0 * y^0)],
    )
    S = typeof(set)
    # `ImageBridge` directly supports `F`-in-`WeightedSOSCone`.
    @test MOI.Bridges.bridging_cost(optimizer, F, S) == 7.0
    func =
        MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[5, -1, 2, 2])
    ci = MOI.add_constraint(optimizer, func, set)
    @test MOI.Bridges.bridge(optimizer, ci) isa
          SumOfSquares.Bridges.Constraint.ImageBridge
    # Re-instantiate without `ImageBridge` to measure the fallback cost. The
    # only remaining route is via a free variable + `KernelBridge`, but the
    # full chain (functionize/scalarize/slack) is strictly more expensive than
    # `ImageBridge`, so `ImageBridge` stays the chosen path when both are
    # available.
    optimizer = MOI.instantiate(
        Clarabel.Optimizer;
        with_cache_type = T,
        with_bridge_type = T,
    )
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    MOI.Bridges.remove_bridge(
        optimizer,
        SumOfSquares.Bridges.Constraint.ImageBridge{T},
    )
    @test MOI.Bridges.bridging_cost(optimizer, F, S) == 8.0
    return
end

end

TestLazy.runtests()
