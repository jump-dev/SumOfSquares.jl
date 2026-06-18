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

end  # module

TestConstraintSOSPolynomialInSemialgebraicSet.runtests()
