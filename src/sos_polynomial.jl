function build_gram_matrix(q::Vector{MOI.VariableIndex},
                           basis::MB.AbstractPolynomialBasis)
    return build_gram_matrix([MOI.SingleVariable(vi) for vi in q], basis)
end
function build_gram_matrix(q::Vector,
                           basis::MB.AbstractPolynomialBasis)
    return GramMatrix(MultivariateMoments.SymMatrix(q, length(basis)),
                      basis)
end

function add_gram_matrix(model::MOI.ModelLike, matrix_cone_type::Type,
                         basis::MB.AbstractPolynomialBasis)
    Q, cQ = MOI.add_constrained_variables(model, matrix_cone(matrix_cone_type, length(basis)))
    q = build_gram_matrix(Q, basis)
    return q, Q, cQ
end
function add_gram_matrix(model::MOI.ModelLike, matrix_cone_type::Type,
                         bases::Vector{<:MB.AbstractPolynomialBasis})
    qQcQ = add_gram_matrix.(model, matrix_cone_type, bases)
    return SparseGramMatrix(getindex.(qQcQ, 1)), getindex.(qQcQ, 2), getindex.(qQcQ, 3)
end

function build_moment_matrix(q::Vector,
                             basis::MB.AbstractPolynomialBasis)
    return MomentMatrix(MultivariateMoments.SymMatrix(q, length(basis)),
                        basis)
end

struct SOSPolynomialSet{DT <: AbstractSemialgebraicSet,
                        MT <: MP.AbstractMonomial,
                        MVT <: AbstractVector{MT},
                        CT <: Certificate.AbstractCertificate} <: MOI.AbstractVectorSet
    domain::DT
    monomials::MVT
    certificate::CT
end
MOI.dimension(set::SOSPolynomialSet) = length(set.monomials)
