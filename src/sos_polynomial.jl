import ComplexOptInterface
const COI = ComplexOptInterface

function matrix_cone(::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, d)
    return COI.HermitianPositiveSemidefiniteConeTriangle(d)
end

vectorized_matrix(Q, basis, ::Type, ::Type) = MultivariateMoments.SymMatrix(Q, basis)
function vectorized_matrix(Q, basis, ::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, T::Type)
    return MultivariateMoments.VectorizedHermitianMatrix{eltype(Q), T}(Q, basis)
end
function vectorized_matrix(Q, basis, ::Type{COI.HermitianPositiveSemidefiniteConeTriangle}, ::Type{Complex{T}}) where T
    return vectorized_matrix(Q, basis, COI.HermitianPositiveSemidefiniteConeTriangle, T)
end

# Need these two methods to avoid ambiguity
function build_gram_matrix(q::Vector{MOI.VariableIndex},
                           basis::AbstractPolynomialBasis, matrix_cone_type, T::Type)
    return build_gram_matrix([MOI.SingleVariable(vi) for vi in q], basis, matrix_cone_type, T)
end
#function build_gram_matrix(q::Vector{MOI.VariableIndex},
#                           basis::AbstractPolynomialBasis, T::Type, vectorization::HermitianVectorized)
#    return build_gram_matrix([MOI.SingleVariable(vi) for vi in q], basis, T, vectorization)
#end
function build_gram_matrix(q::Vector,
                           basis::AbstractPolynomialBasis, matrix_cone_type, T::Type)
    n = length(basis)
    Q = vectorized_matrix(q, length(basis), matrix_cone_type, T)
    U = _promote_sum(eltype(Q), T)
#    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
#    if length(q) != N
#        throw(DimensionMismatch("Length of vectorization $(length(q)) does not match the side dimension $(n) of the symmetric matrix, it should be $N."))
#    end
    return GramMatrix{eltype(Q), typeof(basis), U}(Q, basis)
end
#function build_gram_matrix(q::Vector,
#                           basis::AbstractPolynomialBasis, T::Type, ::HermitianVectorized)
#    n = length(basis)
#    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
#    M = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n - 1))
#    U = _promote_sum(eltype(q), T)
#    if length(q) != N + M
#        throw(DimensionMismatch("Length of real vectorization $(length(q)) does not match the side dimension $(n) of the hermitian matrix, it should be $(N + M) for complex hermitian matrices."))
#    end
#    # Complex case, the first part of `q` is the real part,
#    # the second part is the complex part of the upper off-diagonal entries.
#    C = [zero(U) for i in 1:N]
#    k_real = 0
#    k_imag = 0
#    for j in 1:n
#        for i in 1:j
#            k_real += 1
#            C[k_real] = MA.operate!(MA.add_mul, C[k_real], one(T), q[k_real])
#            if i != j
#                k_imag += 1
#                C[k_real] = MA.operate!(MA.add_mul, C[k_real], one(T) * im, q[N + k_imag])
#            end
#        end
#    end
#    return build_gram_matrix(C, basis, T, SymmetricVectorized())
#end

function union_constraint_indices_types(MCT)
    return Union{MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 0))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 1))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 2))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(matrix_cone(MCT, 3))}}
end

function add_gram_matrix(model::MOI.ModelLike, matrix_cone_type::Type,
                         basis::AbstractPolynomialBasis, T::Type)
    Q, cQ = MOI.add_constrained_variables(model, matrix_cone(matrix_cone_type, length(basis)))
    q = build_gram_matrix(Q, basis, matrix_cone_type, T)
    return q, Q, cQ
end
function add_gram_matrix(model::MOI.ModelLike, matrix_cone_type::Type,
                         bases::Vector{<:AbstractPolynomialBasis}, T::Type)
    Qs = Vector{Vector{MOI.VariableIndex}}(undef, length(bases))
    cQs = Vector{union_constraint_indices_types(matrix_cone_type)}(undef, length(bases))
    # We use `map` for `grams` as it's less easy to infer its type
    grams = map(eachindex(bases)) do i
        gram, Qs[i], cQs[i] = add_gram_matrix(model, matrix_cone_type, bases[i], T)
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
