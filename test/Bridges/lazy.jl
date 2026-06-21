module TestLazy

using Test
import Clarabel
using DynamicPolynomials
import Dualization
using JuMP
import MultivariateBases as MB
import MultivariatePolynomials as MP
using SumOfSquares

# Build a `Putinar(Newton, Newton, maxdegree)` certificate, matching what
# `@constraint(model, ... in SOSCone(), domain = K, maxdegree = d)` produces
# under the hood.
function _newton_putinar(vars, maxdegree)
    full_basis = MB.FullBasis{MB.Monomial}(vars)
    newton = SumOfSquares.Certificate.Newton(
        SOSCone(),
        full_basis,
        full_basis,
        SumOfSquares.Certificate.NewtonFilter(
            SumOfSquares.Certificate.NewtonDegreeBounds(tuple()),
        ),
    )
    return SumOfSquares.Certificate.Putinar(newton, newton, maxdegree)
end

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
    @polyvar x y
    T = Float64
    optimizer = MOI.instantiate(
        Clarabel.Optimizer;
        with_cache_type = T,
        with_bridge_type = T,
    )
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    F = MOI.VectorAffineFunction{T}
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
        MB.SubBasis{MB.Monomial}([y^4, x^2 * y^2, x^3 * y, x^4]),
        [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
        [MB.algebra_element(one(T) * x^0 * y^0)],
    )
    S = typeof(set)
    # `ImageBridge` directly supports `F`-in-`WeightedSOSCone`.
    @test MOI.Bridges.bridging_cost(optimizer, F, S) == 6.0
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

# Wrapping `Clarabel.Optimizer` in `Dualization.dual_optimizer` flips
# Clarabel's `VAF`-in-`PositiveSemidefiniteConeTriangle` support into
# `PositiveSemidefiniteConeTriangle` supported as constrained variables on
# the outer bridge graph. The cost graph propagates the inner solver's
# support through the dualization layer, so the variable-side route wins
# naturally: `VectorSlackBridge` lifts the `VAF`-in-`WeightedSOSCone`
# constraint into constrained variables in `WeightedSOSCone`, which
# `Variable.KernelBridge` then breaks into PSD blocks.
function test_dual_clarabel_uses_kernel_bridge()
    T = Float64
    @polyvar x y
    optimizer = MOI.instantiate(
        Dualization.dual_optimizer(Clarabel.Optimizer);
        with_cache_type = T,
        with_bridge_type = T,
    )
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    set = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
        MB.SubBasis{MB.Monomial}([y^4, x * y^3, x^2 * y^2, x^3 * y, x^4]),
        [MB.SubBasis{MB.Monomial}([y^2, x * y, x^2])],
        [MB.algebra_element(one(T) * x^0 * y^0)],
    )
    func = MOI.VectorAffineFunction{T}(
        MOI.VectorAffineTerm{T}[],
        T[5, -1, 2, 2, 0],
    )
    ci = MOI.add_constraint(optimizer, func, set)
    outer = MOI.Bridges.bridge(optimizer, ci)
    @test outer isa MOI.Bridges.Constraint.VectorSlackBridge
    inner = MOI.Bridges.bridge(optimizer, outer.slack_in_set)
    @test inner isa SumOfSquares.Bridges.Variable.KernelBridge
    return
end

# With a `domain` keyword, the constraint enters the bridge graph as
# `SOSPolynomialSet{<:BasicSemialgebraicSet}`, which
# `SOSPolynomialInSemialgebraicSetBridge` reformulates into a single
# `WeightedSOSCone`. From there, Clarabel — which natively supports
# `MOI.PositiveSemidefiniteConeTriangle` as a constraint — routes through
# `Constraint.ImageBridge`, mirroring `test_clarabel_uses_image_bridge`.
function test_clarabel_with_domain_uses_image_bridge()
    T = Float64
    @polyvar x
    optimizer = MOI.instantiate(
        Clarabel.Optimizer;
        with_cache_type = T,
        with_bridge_type = T,
    )
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    basis = MB.SubBasis{MB.Monomial}([x, x^2])
    domain = (@set x >= 0)
    cert = _newton_putinar([x], 2)
    set = SumOfSquares.SOSPolynomialSet(domain, basis, cert)
    F = MOI.VectorAffineFunction{T}
    S = typeof(set)
    # `SOSPolynomialInSemialgebraicSetBridge` (cost 1) + `ImageBridge` path
    # (cost 6 for Clarabel; same as `test_clarabel_uses_image_bridge`).
    @test MOI.Bridges.bridging_cost(optimizer, F, S) == 7.0
    func = MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[1, 1])
    ci = MOI.add_constraint(optimizer, func, set)
    outer = MOI.Bridges.bridge(optimizer, ci)
    @test outer isa
          SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge
    inner = MOI.Bridges.bridge(optimizer, outer.constraint)
    @test inner isa SumOfSquares.Bridges.Constraint.ImageBridge
    return
end

# Wrapping `Clarabel.Optimizer` in `Dualization.dual_optimizer` flips
# Clarabel's `VAF`-in-`PositiveSemidefiniteConeTriangle` support into
# `PositiveSemidefiniteConeTriangle` supported as constrained variables on
# the outer bridge graph. As in `test_dual_clarabel_uses_kernel_bridge`,
# the cost graph propagates this through the layers and the variable-side
# route wins naturally: `SOSPolynomialInSemialgebraicSetBridge` →
# `VectorSlackBridge` → `Variable.KernelBridge`.
function test_dual_clarabel_with_domain_uses_kernel_bridge()
    T = Float64
    @polyvar x
    optimizer = MOI.instantiate(
        Dualization.dual_optimizer(Clarabel.Optimizer);
        with_cache_type = T,
        with_bridge_type = T,
    )
    SumOfSquares.Bridges.add_all_bridges(optimizer, T)
    basis = MB.SubBasis{MB.Monomial}([x, x^2])
    domain = (@set x >= 0)
    cert = _newton_putinar([x], 2)
    set = SumOfSquares.SOSPolynomialSet(domain, basis, cert)
    func = MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[1, 1])
    ci = MOI.add_constraint(optimizer, func, set)
    outer = MOI.Bridges.bridge(optimizer, ci)
    @test outer isa
          SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge
    slack = MOI.Bridges.bridge(optimizer, outer.constraint)
    @test slack isa MOI.Bridges.Constraint.VectorSlackBridge
    inner = MOI.Bridges.bridge(optimizer, slack.slack_in_set)
    @test inner isa SumOfSquares.Bridges.Variable.KernelBridge
    return
end

# Walk the bridge chain on a JuMP model directly (no `MOI.instantiate`,
# no explicit `add_all_bridges` call). `PolyJuMP.bridges` is responsible
# for transitively registering every bridge the chain might use; if any
# is missing, the lazy cost graph picks the wrong path. The JuMP-side
# chain should match `test_clarabel_with_domain_uses_image_bridge`.
function test_clarabel_with_domain_uses_image_bridge_jump()
    @polyvar x
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, α)
    @objective(model, Max, α)
    S = @set x >= 0
    @constraint(model, α - x^2 in SOSCone(), domain = S, maxdegree = 2)
    optimize!(model)
    lbo = JuMP.backend(model).optimizer
    cis = MOI.ConstraintIndex[]
    for (F, S) in MOI.get(lbo, MOI.ListOfConstraintTypesPresent())
        if S <: SumOfSquares.SOSPolynomialSet
            append!(cis, MOI.get(lbo, MOI.ListOfConstraintIndices{F,S}()))
        end
    end
    @test length(cis) == 1
    outer = MOI.Bridges.bridge(lbo, only(cis))
    @test outer isa
          SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge
    inner = MOI.Bridges.bridge(lbo, outer.constraint)
    @test inner isa SumOfSquares.Bridges.Constraint.ImageBridge
    return
end

# JuMP counterpart of `test_dual_clarabel_with_domain_uses_kernel_bridge`:
# wrap Clarabel in `Dualization.dual_optimizer` and verify the chain is
# `SOSPolynomialInSemialgebraicSetBridge` → `VectorSlackBridge` →
# `Variable.KernelBridge` without any explicit bridge removal.
function test_dual_clarabel_with_domain_uses_kernel_bridge_jump()
    @polyvar x
    model = Model(Dualization.dual_optimizer(Clarabel.Optimizer))
    set_silent(model)
    @variable(model, α)
    @objective(model, Max, α)
    S = @set x >= 0
    @constraint(model, α - x^2 in SOSCone(), domain = S, maxdegree = 2)
    optimize!(model)
    lbo = JuMP.backend(model).optimizer
    cis = MOI.ConstraintIndex[]
    for (F, S) in MOI.get(lbo, MOI.ListOfConstraintTypesPresent())
        if S <: SumOfSquares.SOSPolynomialSet
            append!(cis, MOI.get(lbo, MOI.ListOfConstraintIndices{F,S}()))
        end
    end
    @test length(cis) == 1
    outer = MOI.Bridges.bridge(lbo, only(cis))
    @test outer isa
          SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge
    slack = MOI.Bridges.bridge(lbo, outer.constraint)
    @test slack isa MOI.Bridges.Constraint.VectorSlackBridge
    inner = MOI.Bridges.bridge(lbo, slack.slack_in_set)
    @test inner isa SumOfSquares.Bridges.Variable.KernelBridge
    return
end

end

TestLazy.runtests()
