function build_gram_matrix(q::Vector,
                           monos::AbstractVector{<:MP.AbstractMonomial})
    return GramMatrix(MultivariateMoments.SymMatrix(q, length(monos)),
                         monos)
end

function build_moment_matrix(q::Vector,
                             monos::AbstractVector{<:MP.AbstractMonomial})
    return MomentMatrix(MultivariateMoments.SymMatrix(q, length(monos)),
                        monos)
end

abstract type SOSLikeCone <: PolyJuMP.PolynomialSet end
Base.broadcastable(cone::SOSLikeCone) = Ref(cone)

struct SOSPolynomialSet{DT <: AbstractSemialgebraicSet,
                        CT <: SOSLikeCone,
                        BT <: PolyJuMP.AbstractPolynomialBasis,
                        MT <: MP.AbstractMonomial,
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
