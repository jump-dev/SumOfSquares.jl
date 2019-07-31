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

struct SOSPolynomialSet{DT <: AbstractSemialgebraicSet,
                        MT <: MP.AbstractMonomial,
                        MVT <: AbstractVector{MT},
                        CT <: Certificate.AbstractCertificate} <: MOI.AbstractVectorSet
    domain::DT
    monomials::MVT
    certificate::CT
end
