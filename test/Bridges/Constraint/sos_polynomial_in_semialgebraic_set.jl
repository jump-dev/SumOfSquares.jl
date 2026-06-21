module TestConstraintSOSPolynomialInSemialgebraicSet

using Test
import MathOptInterface as MOI
import MultivariateBases as MB
import MultivariatePolynomials as MP
import StarAlgebras as SA
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

# Build a `Putinar(Newton, Newton, maxdegree)` certificate, the default that
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

# Build a `Sparsity.Preorder(Variable, Putinar(Newton, MaxDegree, maxdegree))`
# certificate, matching what `sparsity = Sparsity.Variable()` produces in
# JuMP.
function _sparse_putinar(vars, maxdegree)
    full_basis = MB.FullBasis{MB.Monomial}(vars)
    newton = SumOfSquares.Certificate.Newton(
        SOSCone(),
        full_basis,
        full_basis,
        SumOfSquares.Certificate.NewtonFilter(
            SumOfSquares.Certificate.NewtonDegreeBounds(tuple()),
        ),
    )
    max_deg = SumOfSquares.Certificate.MaxDegree(
        SOSCone(),
        full_basis,
        full_basis,
        maxdegree,
    )
    return SumOfSquares.Certificate.Sparsity.Preorder(
        SumOfSquares.Certificate.Sparsity.Variable(),
        SumOfSquares.Certificate.Putinar(newton, max_deg, maxdegree),
    )
end

# Cover the core code path of the bridge: σ_0's gram basis is determined by
# `Newton` on the *augmented* polynomial `p + Σ g_i · b_i ⋆ b_i` (a recent
# bugfix; without the augmentation, Newton on `p` alone gives a different
# basis), and a single Lagrangian multiplier σ_1 is added for the
# inequality `x ≥ 0`.
#
# Input polynomial `p(x) = (1 - a) x + a x^2` with `basis = [x, x^2]`,
# `domain = {x ≥ 0}` and `maxdegree = 2`.
# - multiplier basis for `g_1 = x`: `[1]`
# - augmented poly monomials: `{x, x^2}` (Newton polytope ⇒ gram basis `[x]`)
# - σ_0 gram basis: `[x]`
# - σ_1 gram basis: `[1]`
# - WeightedSOSCone weights: `[g_1, 1] = [x, 1]`
# - `new_basis = [x, x^2]` (coefficients align with the input function).
function test_runtests()
    T = Float64
    @polyvar x
    basis = MB.SubBasis{MB.Monomial}([x, x^2])
    domain = (@set x >= 0)
    cert = _newton_putinar([x], 2)
    set = SumOfSquares.SOSPolynomialSet(domain, basis, cert)
    sigma_0_basis = MB.SubBasis{MB.Monomial}([x])
    sigma_1_basis = MB.SubBasis{MB.Monomial}([x^0])
    new_basis = MB.SubBasis{MB.Monomial}([x, x^2])
    implicit_basis = MB.implicit_basis(basis)
    w_unit = MB.algebra_element(
        MB.sparse_coefficients(
            MP.polynomial(MP.term(one(T), MP.constant_monomial(x))),
        ),
        implicit_basis,
    )
    g1 = MP.polynomial(domain.p[1])
    w_g1 = MB.algebra_element(
        MB.sparse_coefficients(one(T) * MP.similar(g1, T)),
        implicit_basis,
    )
    weighted =
        SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
            new_basis,
            [sigma_1_basis, sigma_0_basis],
            [w_g1, w_unit],
        )
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge,
        model -> begin
            a = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.Utilities.operate(
                    vcat,
                    T,
                    T(1) - T(1) * a,
                    T(0) + T(1) * a,
                ),
                set,
            )
        end,
        model -> begin
            a = MOI.add_variable(model)
            MOI.add_constraint(
                model,
                MOI.Utilities.operate(
                    vcat,
                    T,
                    T(1) - T(1) * a,
                    T(0) + T(1) * a,
                ),
                weighted,
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# Two inequalities to exercise the `multiplier_indices` machinery for σ_0
# (last) and σ_1, σ_2 (first two `gram_bases`). The augmentation now adds
# monomials from *both* `g_i b_i ⋆ b_i` products to the polynomial used to
# compute σ_0's Newton polytope.
function test_runtests_two_inequalities()
    T = Float64
    @polyvar x y
    basis = MB.SubBasis{MB.Monomial}([x * y])
    domain = @set x >= 0 && y >= 0
    cert = _newton_putinar([x, y], 2)
    set = SumOfSquares.SOSPolynomialSet(domain, basis, cert)
    # With the constant input function `1`, Newton's filter on the augmented
    # polynomial gives empty multiplier gram bases — what matters for this
    # test is that we exercise the multi-inequality `for index in preorder_indices`
    # loop and the `multiplier_indices` bookkeeping. We have to use the
    # same parent `FullBasis{Monomial}([x, y])` as the certificate's gram
    # basis (the parent is part of the `SubBasis`'s identity).
    sigma_basis =
        SA.SubBasis(MB.FullBasis{MB.Monomial}([x, y]), Vector{Vector{Int}}())
    new_basis = MB.SubBasis{MB.Monomial}([x * y])
    implicit_basis = MB.implicit_basis(basis)
    w_unit = MB.algebra_element(
        MB.sparse_coefficients(
            MP.polynomial(MP.term(one(T), MP.constant_monomial(x * y))),
        ),
        implicit_basis,
    )
    w_g1 = MB.algebra_element(
        MB.sparse_coefficients(
            one(T) * MP.similar(MP.polynomial(domain.p[1]), T),
        ),
        implicit_basis,
    )
    w_g2 = MB.algebra_element(
        MB.sparse_coefficients(
            one(T) * MP.similar(MP.polynomial(domain.p[2]), T),
        ),
        implicit_basis,
    )
    weighted =
        SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
            new_basis,
            [sigma_basis, sigma_basis, sigma_basis],
            [w_g1, w_g2, w_unit],
        )
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[1]),
                set,
            )
        end,
        model -> begin
            MOI.add_constraint(
                model,
                MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[1]),
                weighted,
            )
        end;
        cannot_unbridge = true,
    )
    return
end

# Exercise the `if raw isa AbstractVector` branch of the `raw_ranges`
# bookkeeping: when the certificate carries sparsity (here
# `Sparsity.Preorder(Variable, Putinar(Newton, MaxDegree, …))`), both the σ_i
# multiplier bases and σ_0's basis come back as `Vector{SubBasis}` and the
# bridge collapses them into a `UnitRange` slot of `multiplier_indices`. For
# `p = x^2 + y^2` over `domain = {x ≥ 0}`, the σ_0 basis is split into two
# sparsity blocks `[1, x]` and `[1, y]` (one per disconnected variable
# cluster), so `multiplier_indices = [2:3, 1]`.
function test_runtests_sparsity()
    T = Float64
    @polyvar x y
    basis = MB.SubBasis{MB.Monomial}([y^2, x^2])
    domain = (@set x >= 0)
    cert = _sparse_putinar([x, y], 2)
    set = SumOfSquares.SOSPolynomialSet(domain, basis, cert)
    full = MB.FullBasis{MB.Monomial}([x, y])
    # Construct the basis and gram bases with the same parent
    # `FullBasis{Monomial}([x, y])` as the certificate uses, so that
    # `==` on `SubBasis` returns true.
    new_basis = SA.SubBasis(full, [[0, 0], [0, 1], [1, 0], [0, 2], [2, 0]])
    sigma_1_basis = SA.SubBasis(full, [[0, 0]])
    sigma_0_block_y = SA.SubBasis(full, [[0, 0], [0, 1]])
    sigma_0_block_x = SA.SubBasis(full, [[0, 0], [1, 0]])
    implicit_basis = MB.implicit_basis(basis)
    w_unit = MB.algebra_element(
        MB.sparse_coefficients(
            MP.polynomial(MP.term(one(T), MP.constant_monomial(x * y))),
        ),
        implicit_basis,
    )
    # `domain.p[1] = x` only carries the variable `x`, so we lift it into
    # the `[x, y]` variable system with `+ 0*y` to match the bridge's
    # `Certificate.generator(...)` output exactly.
    g1_full = MP.polynomial(domain.p[1] + 0 * y)
    w_g1 = MB.algebra_element(
        MB.sparse_coefficients(one(T) * MP.similar(g1_full, T)),
        implicit_basis,
    )
    weighted =
        SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
            new_basis,
            [sigma_1_basis, sigma_0_block_x, sigma_0_block_y],
            [w_g1, w_unit, w_unit],
        )
    func = MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[1, 1])
    new_func =
        MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[0, 0, 0, 1, 1])
    MOI.Bridges.runtests(
        SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge,
        model -> begin
            MOI.add_constraint(model, func, set)
        end,
        model -> begin
            MOI.add_constraint(model, new_func, weighted)
        end;
        cannot_unbridge = true,
    )
    return
end

# Verify the `multiplier_indices` shape that the `GramMatrixAttribute`,
# `MomentMatrixAttribute` and `LagrangianMultipliers` getters dispatch on:
# - σ_0 at `multiplier_indices[1]` must be a `UnitRange` (the two sparsity
#   blocks ⇒ `_get` ⇒ `BlockDiagonalGramMatrix`);
# - σ_1 at `multiplier_indices[2]` must also be a `UnitRange` (single
#   block wrapped in a `Vector` ⇒ same `BlockDiagonalGramMatrix` path).
# The actual gram-matrix recovery on a real-solver path is covered by
# `test/Tests/rearrangement.jl` (run from `test/Mock/rearrangement.jl`).
function test_multiplier_indices_sparsity()
    T = Float64
    @polyvar x y
    basis = MB.SubBasis{MB.Monomial}([y^2, x^2])
    domain = (@set x >= 0)
    cert = _sparse_putinar([x, y], 2)
    set = SumOfSquares.SOSPolynomialSet(domain, basis, cert)
    # Build the bridged model with a mock optimizer that we can poke.
    mock = MOI.Utilities.MockOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{T}()),
    )
    wrapped = MOI.Bridges.Constraint.SingleBridgeOptimizer{
        SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge{
            T,
        },
    }(
        mock,
    )
    func = MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], T[1, 1])
    ci = MOI.add_constraint(wrapped, func, set)
    bridge = MOI.Bridges.bridge(wrapped, ci)
    # σ_0 occupies `multiplier_indices[1]` and is split into 2 sparsity
    # blocks, hence the `UnitRange` `2:3`. σ_1 has a single sparsity block
    # but the sparsity certificate still wraps its basis in a `Vector`, so
    # it comes back as the (length-1) `UnitRange` `1:1` instead of a plain
    # `Int` — that single-block-but-Vector path is precisely the case the
    # `BlockDiagonalGramMatrix` getter has to handle.
    @test bridge.multiplier_indices[1] == 2:3
    @test bridge.multiplier_indices[2] == 1:1
    return
end

end  # module

TestConstraintSOSPolynomialInSemialgebraicSet.runtests()
