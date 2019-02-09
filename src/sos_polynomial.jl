export DSOSCone, SDSOSCone, SOSCone

export CoDSOSCone, CoSDSOSCone, CoSOSCone
export SOSMatrixCone

struct DSOSCone <: PolyJuMP.PolynomialSet end
matrix_cone(::DSOSCone) = DiagonallyDominantConeTriangle

struct SDSOSCone <: PolyJuMP.PolynomialSet end
matrix_cone(::SDSOSCone) = ScaledDiagonallyDominantConeTriangle

struct SOSCone <: PolyJuMP.PolynomialSet end
matrix_cone(::SOSCone) = MOI.PositiveSemidefiniteConeTriangle

const SOSLikeCones = Union{DSOSCone, SDSOSCone, SOSCone}

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
and the cone of symmetric matrices with nonnegative entries
[Lemma 3.164, BPT12].

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
