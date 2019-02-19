export DSOSCone, SDSOSCone, SOSCone
export CopositiveInner
export SOSMatrixCone

function build_gram_matrix(q::Vector,
                           monos::AbstractVector{<:AbstractMonomial})
    return MatPolynomial(MultivariateMoments.SymMatrix(q, length(monos)),
                         monos)
end

function build_moment_matrix(q::Vector,
                             monos::AbstractVector{<:AbstractMonomial})
    return MomentMatrix(MultivariateMoments.SymMatrix(q, length(monos)),
                        monos)
end

abstract type SOSLikeCone <: PolyJuMP.PolynomialSet end
Base.broadcastable(cone::SOSLikeCone) = Ref(cone)

abstract type SOSSubCone <: SOSLikeCone end

struct DSOSCone <: SOSSubCone end
matrix_cone_type(::Type{DSOSCone}) = DiagonallyDominantConeTriangle

struct SDSOSCone <: SOSSubCone end
matrix_cone_type(::Type{SDSOSCone}) = ScaledDiagonallyDominantConeTriangle

struct SOSCone <: SOSSubCone end
matrix_cone_type(::Type{SOSCone}) = MOI.PositiveSemidefiniteConeTriangle

struct SOSPolynomialSet{DT <: AbstractSemialgebraicSet,
                        CT <: SOSLikeCone,
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
