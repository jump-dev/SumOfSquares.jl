function build_gram_matrix(q::Vector{MOI.VariableIndex},
                           basis::AbstractPolynomialBasis)
    return build_gram_matrix([MOI.SingleVariable(vi) for vi in q], basis)
end
function build_gram_matrix(q::Vector,
                           basis::AbstractPolynomialBasis)
    return GramMatrix(MultivariateMoments.SymMatrix(q, length(basis)),
                      basis)
end

function union_constraint_indices_types(MCT)
    return Union{MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 0))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 1))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 2))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 3))}}
end

function add_gram_matrix(model::MOI.ModelLike, matrix_cone_type::Type,
                         basis::AbstractPolynomialBasis)
    Q, cQ = MOI.add_constrained_variables(model, matrix_cone(matrix_cone_type, length(basis)))
    q = build_gram_matrix(Q, basis)
    return q, Q, cQ
end
function add_gram_matrix(model::MOI.ModelLike, matrix_cone_type::Type,
                         bases::Vector{<:AbstractPolynomialBasis})
    Qs = Vector{Vector{MOI.VariableIndex}}(undef, length(bases))
    cQs = Vector{union_constraint_indices_types(matrix_cone_type)}(undef, length(bases))
    # We use `map` for `grams` as it's less easy to infer its type
    grams = map(eachindex(bases)) do i
        gram, Qs[i], cQs[i] = add_gram_matrix(model, matrix_cone_type, bases[i])
        return gram
    end
    return SparseGramMatrix(grams), Qs, cQs
end

function build_moment_matrix(q::Vector,
                             basis::AbstractPolynomialBasis)
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
Base.copy(set::SOSPolynomialSet) = set
