export DSOSCone, SDSOSCone, SOSCone
export CopositiveInner
export SOSMatrixCone

function add_gram_matrix(model::MOI.ModelLike, monos)
    return MatPolynomial{MOI.SingleVariable}(
        (i, j) -> MOI.SingleVariable(MOI.add_variable(model)), monos)
end
function gram_delete(model::MOI.ModelLike, p::MatPolynomial)
    for sv in p.Q.Q
        MOI.delete(model, sv.variable)
    end
end

struct DSOSCone <: PolyJuMP.PolynomialSet end
matrix_cone(::DSOSCone) = DiagonallyDominantConeTriangle

struct SDSOSCone <: PolyJuMP.PolynomialSet end
matrix_cone(::SDSOSCone) = ScaledDiagonallyDominantConeTriangle

struct SOSCone <: PolyJuMP.PolynomialSet end
matrix_cone(::SOSCone) = MOI.PositiveSemidefiniteConeTriangle

const SOSLikeCones = Union{DSOSCone, SDSOSCone, SOSCone}

function gram_in_cone(model::MOI.ModelLike, monos, set::SOSLikeCones)
    p = add_gram_matrix(model, monos)
    ci = matrix_add_constraint(model, p, matrix_cone(set))
    return p, ci
end

"""
    struct CopositiveInner{S} <: PolyJuMP.PolynomialSet
        # Inner approximation of the PSD cone, i.e. typically either
        # `SOSCone`, `DSOSCone` or `SDSOSCone`,
        psd_inner::S
    end

A symmetric matrix ``Q`` is copositive if ``x^{\\top} Q x \\ge 0`` for all
vector ``x`` in the nonnegative orthant. Checking copositivity is a
co-NP-complete problem [MK87] and this cone is only the inner approximation of
the cone of copositive symmetric matrices given by Minknowski sum of `psd_inner`
and the cone of symmetric matrices with nonnegative entries (the diagonal
entries can be chosen to be zero) [Lemma 3.164, BPT12].

The matrix with nonnegative entries can be interpreted as lagrangian
multipliers. For instance,
```julia
@polyvar x y
@constraint(model, x^2 - 2x*y + y^2 in CopositiveInner(SOSCone()))
```
is equivalent to
```julia
# Matrix that we require to be copositive
Q = [ 1 -1
     -1  1]
λ = @variable(model, lower_bound=0)
# Symmetric matrix of nonnegative entries
Λ = [0 λ
     λ 0]
using LinearAlgebra # For `Symmetric`
@constraint(model, Symmetric(Q - Λ) in PSDCone())
```
which is equivalent to
```julia
@polyvar x y
λ = @variable(model, lower_bound=0)
@constraint(model, x^2 - 2x*y + y^2 - 2*λ * x*y in SOSCone())
```
which is the same as, using the `domain` keyword,
```julia
@polyvar x y
@constraint(model, x^2 - 2x*y + y^2 in SOSCone(), domain = @set x*y ≥ 0)
```

For consistency with its equivalent forms, the [`GramMonomials`](@ref) for this
constraint is given by the gram matrix in the `psd_inner` cone, i.e. which
should be equal to `Q - Λ`.

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.

[MK87] K. G. Murty and S. N. Kabadi.
*Some NP-complete problems in quadratic and nonlinear programming*.
Mathematical programming, 39:117–129, **1987**.
"""
struct CopositiveInner{S} <: PolyJuMP.PolynomialSet
    psd_inner::S
end

function gram_posynomial(model::MOI.ModelLike, monos)
    # TODO, the diagonal elements can be zero
    p = add_gram_matrix(model, monos)
    # TODO use Nonnegatives cone
    for q in p.Q.Q
        MOI.add_constraint(model, q, MOI.GreaterThan(0.0))
    end
    return p
end

function gram_in_cone(model::MOI.ModelLike, x, set::CopositiveInner)
    p, ci = gram_in_cone(model, x, set.psd_inner)
    _matplus(p, gram_posynomial(model, x)), ci
end

const SOSSubCones = Union{CopositiveInner, SOSLikeCones}

Base.broadcastable(cone::SOSSubCones) = Ref(cone)

struct SOSPolynomialSet{DT <: AbstractSemialgebraicSet,
                        CT <: SOSSubCones,
                        BT <: PolyJuMP.AbstractPolynomialBasis,
                        MT <: AbstractMonomial,
                        MVT <: AbstractVector{MT},
                        NPT <: Tuple} <: MOI.AbstractVectorSet
    domain::DT
    cone::CT
    basis::Type{BT}
    monomials::MVT
    newton_polytope::NPT
    mindegree::Int
    maxdegree::Int
end
