export GramMatrix, BlockDiagonalGramMatrix

import MultivariateMoments: vectorized_symmetric_matrix, SymMatrix, value_matrix
export gram_operate, value_matrix

abstract type AbstractDecomposition{T} <: MP.AbstractPolynomialLike{T} end

Base.:(+)(x::_APL, y::AbstractDecomposition) = x + MP.polynomial(y)
Base.:(+)(x::AbstractDecomposition, y::_APL) = MP.polynomial(x) + y
function Base.:(+)(x::AbstractDecomposition, y::AbstractDecomposition)
    return MP.polynomial(x) + MP.polynomial(y)
end
Base.:(-)(x::_APL, y::AbstractDecomposition) = x - MP.polynomial(y)
Base.:(-)(x::AbstractDecomposition, y::_APL) = MP.polynomial(x) - y
function Base.:(-)(x::AbstractDecomposition, y::AbstractDecomposition)
    return MP.polynomial(x) - MP.polynomial(y)
end

abstract type AbstractGramMatrix{T,B,U} <: AbstractDecomposition{U} end

function MP.monomial_type(
    ::Union{AbstractGramMatrix{T,B},Type{<:AbstractGramMatrix{T,B}}},
) where {T,B}
    return MP.monomial_type(B)
end
function MP.term_type(
    p::Union{AbstractGramMatrix{T,B,U},Type{<:AbstractGramMatrix{T,B,U}}},
) where {T,B,U}
    return MP.term_type(MP.polynomial_type(p))
end
function MP.polynomial_type(
    ::Union{AbstractGramMatrix{T,B,U},Type{<:AbstractGramMatrix{T,B,U}}},
) where {T,B,U}
    return MP.polynomial_type(B, U)
end

Base.:(==)(p::_APL, q::AbstractGramMatrix) = p == MP.polynomial(q)
Base.:(==)(p::AbstractGramMatrix, q::AbstractGramMatrix) = iszero(p - q)

"""
    struct GramMatrix{T, B, U, MT <: AbstractMatrix{T}} <: AbstractGramMatrix{T, B, U}
        Q::MT
        basis::B
    end

Gram matrix ``x^\\top Q x`` where `Q` is a symmetric matrix indexed by the
vector of polynomials of the basis `basis`.
"""
struct GramMatrix{T,B,U,MT<:AbstractMatrix{T}} <: AbstractGramMatrix{T,B,U}
    Q::MT
    basis::B
end
function GramMatrix{T,B,U}(
    Q::AbstractMatrix{T},
    basis::B,
) where {T,B<:AbstractPolynomialBasis,U}
    return GramMatrix{T,B,U,typeof(Q)}(Q, basis)
end
function GramMatrix{T,B}(
    Q::AbstractMatrix{T},
    basis::B,
) where {T,B<:AbstractPolynomialBasis}
    return GramMatrix{T,B,_promote_sum(T)}(Q, basis)
end
function GramMatrix(
    # We don't use `AbstractMatrix` to avoid clash with method below with `σ`.
    Q::Union{
        MultivariateMoments.SymMatrix{T},
        MultivariateMoments.VectorizedHermitianMatrix{T},
    },
    basis::AbstractPolynomialBasis,
) where {T}
    return GramMatrix{T,typeof(basis)}(Q, basis)
end

function MP.similar_type(
    ::Type{GramMatrix{T,B,U,MT}},
    ::Type{S},
) where {T,B,U,MT,S}
    US = _promote_sum(S)
    MS = MP.similar_type(MT, S)
    return GramMatrix{S,B,US,MS}
end

# When taking the promotion of a GramMatrix of JuMP.Variable with a Polynomial JuMP.Variable, it should be a Polynomial of AffExpr
MP.constant_monomial(p::GramMatrix) = MP.constant_monomial(MP.monomial_type(p))
MP.variables(p::GramMatrix) = MP.variables(p.basis)
MP.nvariables(p::GramMatrix) = MP.nvariables(p.basis)

function Base.zero(::Type{GramMatrix{T,B,U,MT}}) where {T,B,U,MT}
    return GramMatrix{T,B,U,MT}(
        similar(MT, (0, 0)),
        MultivariateBases.empty_basis(B),
    )
end
Base.iszero(p::GramMatrix) = iszero(MP.polynomial(p))

Base.getindex(p::GramMatrix, I...) = getindex(p.Q, I...)
Base.copy(p::GramMatrix) = GramMatrix(copy(p.Q), copy(p.basis))

MultivariateMoments.value_matrix(p::GramMatrix{T}) where {T} = p.Q

function GramMatrix{T}(
    f::Function,
    basis::AbstractPolynomialBasis,
    σ = 1:length(basis),
) where {T}
    return GramMatrix{T,typeof(basis)}(
        vectorized_symmetric_matrix(T, f, length(basis), σ),
        basis,
    )
end
function GramMatrix{T}(f::Function, monos::AbstractVector) where {T}
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return GramMatrix{T}(f, MonomialBasis(sorted_monos), σ)
end

function GramMatrix(
    Q::AbstractMatrix{T},
    basis::AbstractPolynomialBasis,
    σ = 1:length(basis),
) where {T}
    return GramMatrix{T}((i, j) -> Q[σ[i], σ[j]], basis)
end
function GramMatrix(Q::AbstractMatrix, monos::AbstractVector)
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return GramMatrix(Q, MonomialBasis(sorted_monos), σ)
end

#function Base.convert{T, PT <: AbstractPolynomial{T}}(::Type{PT}, p::GramMatrix)
#    # coefficient_type(p) may be different than T and MP.polynomial(p) may be different than PT (different module)
#    convert(PT, MP.polynomial(p))
#end
function MP.polynomial(p::GramMatrix{T,B,U}) where {T,B,U}
    return MP.polynomial(p, U)
end
function MP.polynomial(p::GramMatrix, ::Type{S}) where {S}
    return MP.polynomial(value_matrix(p), p.basis, S)
end

function change_basis(
    p::GramMatrix{T,B},
    ::Type{B},
) where {T,B<:AbstractPolynomialBasis}
    return p
end
function change_basis(p::GramMatrix, B::Type{<:AbstractPolynomialBasis})
    return GramMatrix(MultivariateBases.change_basis(p.Q, p.basis, B)...)
end

"""
    gram_operate(::typeof(+), p::GramMatrix, q::GramMatrix)

Computes the Gram matrix equal to the sum between `p` and `q`. On the opposite,
`p + q` gives a polynomial equal to `p + q`. The polynomial `p + q` can also be
obtained by `polynomial(gram_operate(+, p, q))`.
"""
function gram_operate(
    ::typeof(+),
    p::GramMatrix{S,B,US,SymMatrix{S}},
    q::GramMatrix{T,B,UT,SymMatrix{T}},
) where {S,US,T,UT,B}
    basis, Ip, Iq = MultivariateBases.merge_bases(p.basis, q.basis)
    U = MA.promote_operation(+, S, T)
    n = length(basis)
    Qvec = Vector{U}(undef, div(n * (n + 1), 2))
    MA.operate!(zero, Qvec)
    Q = SymMatrix(Qvec, n)
    for j in 1:n
        for i in 1:j
            if !iszero(Ip[j]) && !iszero(Ip[i])
                MultivariateMoments.symmetric_setindex!(
                    Q,
                    MA.add!!(Q[i, j], p.Q[Ip[i], Ip[j]]),
                    i,
                    j,
                )
            end
            if !iszero(Iq[j]) && !iszero(Iq[i])
                MultivariateMoments.symmetric_setindex!(
                    Q,
                    MA.add!!(Q[i, j], q.Q[Iq[i], Iq[j]]),
                    i,
                    j,
                )
            end
        end
    end
    return GramMatrix(Q, basis)
end
function gram_operate(
    ::typeof(+),
    p::GramMatrix{S,Bp},
    q::GramMatrix{T,Bq},
) where {S,T,Bp,Bq}
    B = promote_type(Bp, Bq)
    return gram_operate(+, change_basis(p, B), change_basis(q, B))
end

"""
    gram_operate(/, p::GramMatrix, α)

Computes the Gram matrix equal to `p / α`. On the opposite,
`p / α` gives a polynomial equal to `p / α`. The polynomial `p / α` can also be
obtained by `polynomial(gram_operate(/, p, α))`.
"""
function gram_operate(::typeof(/), q::GramMatrix, α)
    return GramMatrix(map(x -> x / α, q.Q), q.basis)
end

function Base.isapprox(p::GramMatrix, q::GramMatrix; kwargs...)
    return p.basis == q.basis && isapprox(p.Q, q.Q; kwargs...)
end

function Base.show(io::IO, M::GramMatrix)
    print(io, "GramMatrix")
    show_basis_indexed_matrix(io, M)
end

struct BlockDiagonalGramMatrix{T,B,U,MT} <: AbstractGramMatrix{T,B,U}
    blocks::Vector{GramMatrix{T,B,U,MT}}
end

function _sparse_type(::Type{GramMatrix{T,B,U,MT}}) where {T,B,U,MT}
    return BlockDiagonalGramMatrix{T,B,U,MT}
end

function MP.similar_type(
    ::Type{BlockDiagonalGramMatrix{T,B,U,MT}},
    ::Type{S},
) where {T,B,U,MT,S}
    return _sparse_type(MP.similar_type(GramMatrix{T,B,U,MT}, S))
end

function Base.zero(::Type{BlockDiagonalGramMatrix{T,B,U,MT}}) where {T,B,U,MT}
    return BlockDiagonalGramMatrix(GramMatrix{T,B,U,MT}[])
end
function MP.polynomial(p::BlockDiagonalGramMatrix)
    return mapreduce(
        identity,
        MA.add!!,
        p.blocks,
        init = zero(MP.polynomial_type(p)),
    )
end

function Base.show(io::IO, M::BlockDiagonalGramMatrix)
    print(io, "BlockDiagonalGramMatrix")
    MultivariateMoments.show_basis_indexed_blocks(io, M.blocks)
    return
end
