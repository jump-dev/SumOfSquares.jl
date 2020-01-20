export GramMatrix, SparseGramMatrix

import MultivariateMoments: trimat, SymMatrix, getmat
export gram_operate, getmat

abstract type AbstractGramMatrix{T, B, U} <: MP.AbstractPolynomialLike{U} end

function MP.monomialtype(::Union{AbstractGramMatrix{T, B},
                                 Type{<:AbstractGramMatrix{T, B}}}) where {T, B}
    return MP.monomialtype(B)
end
MP.polynomialtype(::Union{AbstractGramMatrix{T, B, U}, Type{<:AbstractGramMatrix{T, B, U}}}) where {T, B, U} = MP.polynomialtype(B, U)

Base.:(+)(x::MP.APL, y::AbstractGramMatrix) = x + MP.polynomial(y)
Base.:(+)(x::AbstractGramMatrix, y::MP.APL) = MP.polynomial(x) + y
Base.:(+)(x::AbstractGramMatrix, y::AbstractGramMatrix) = MP.polynomial(x) + MP.polynomial(y)
Base.:(-)(x::MP.APL, y::AbstractGramMatrix) = x - MP.polynomial(y)
Base.:(-)(x::AbstractGramMatrix, y::MP.APL) = MP.polynomial(x) - y
Base.:(-)(x::AbstractGramMatrix, y::AbstractGramMatrix) = MP.polynomial(x) - MP.polynomial(y)
Base.:(==)(p::MP.APL, q::AbstractGramMatrix) = p == MP.polynomial(q)
Base.:(==)(p::AbstractGramMatrix, q::AbstractGramMatrix) = iszero(p - q)

"""
    struct GramMatrix{T, B} <: MP.APL{T}
        Q::SymMatrix{T}
        basis::B
    end

Gram matrix ``x^\\top Q x`` where `Q` is a symmetric matrix indexed by the
vector of polynomials of the basis `basis`.
"""
struct GramMatrix{T, B, U} <: AbstractGramMatrix{T, B, U}
    Q::SymMatrix{T}
    basis::B
end
GramMatrix{T, B}(Q::SymMatrix{T}, basis::B) where {T, B<:MB.AbstractPolynomialBasis} = GramMatrix{T, B, _promote_sum(T)}(Q, basis)
function GramMatrix(Q::SymMatrix{T}, basis::MB.AbstractPolynomialBasis) where T
    return GramMatrix{T, typeof(basis)}(Q, basis)
end
# When taking the promotion of a GramMatrix of JuMP.Variable with a Polynomial JuMP.Variable, it should be a Polynomial of AffExpr
MP.constantmonomial(p::GramMatrix) = MP.constantmonomial(MP.monomialtype(p))
MP.variables(p::GramMatrix) = MP.variables(p.basis)
MP.nvariables(p::GramMatrix) = MP.nvariables(p.basis)

Base.zero(::Type{GramMatrix{T, B, U}}) where {T, B, U} = GramMatrix{T, B, U}(SymMatrix{T}(T[], 0), MB.empty_basis(B))
Base.iszero(p::GramMatrix) = iszero(MP.polynomial(p))

Base.getindex(p::GramMatrix, I...) = getindex(p.Q, I...)
Base.copy(p::GramMatrix) = GramMatrix(copy(p.Q), copy(p.basis))

MultivariateMoments.getmat(p::GramMatrix{T}) where {T} = p.Q

function GramMatrix{T}(f::Function, basis::MB.AbstractPolynomialBasis, σ=1:length(basis)) where T
    GramMatrix{T, typeof(basis)}(trimat(T, f, length(basis), σ), basis)
end
function GramMatrix{T}(f::Function, monos::AbstractVector) where T
    σ, sorted_monos = MP.sortmonovec(monos)
    return GramMatrix{T}(f, MB.MonomialBasis(sorted_monos), σ)
end

function GramMatrix(Q::AbstractMatrix{T}, basis::MB.AbstractPolynomialBasis, σ) where T
    return GramMatrix{T}((i, j) -> Q[σ[i], σ[j]], basis)
end
function GramMatrix(Q::AbstractMatrix, monos::AbstractVector)
    σ, sorted_monos = MP.sortmonovec(monos)
    return GramMatrix(Q, MB.MonomialBasis(sorted_monos), σ)
end

#function Base.convert{T, PT <: AbstractPolynomial{T}}(::Type{PT}, p::GramMatrix)
#    # coefficienttype(p) may be different than T and MP.polynomial(p) may be different than PT (different module)
#    convert(PT, MP.polynomial(p))
#end
function MP.polynomial(p::GramMatrix{T, B, U}) where {T, B, U}
    return MP.polynomial(p, U)
end
function MP.polynomial(p::GramMatrix, ::Type{S}) where S
    return MP.polynomial(getmat(p), p.basis, S)
end

"""
    gram_operate(::typeof(+), p::GramMatrix, q::GramMatrix)

Computes the Gram matrix equal to the sum between `p` and `q`. On the opposite,
`p + q` gives a polynomial equal to `p + q`. The polynomial `p + q` can also be
obtained by `polynomial(gram_sum(p, q))`.
"""
function gram_operate(::typeof(+), p::GramMatrix{S}, q::GramMatrix{T}) where {S, T}
    basis, Ip, Iq = MB.merge_bases(p.basis, q.basis)
    U = MA.promote_operation(+, S, T)
    n = length(basis)
    Qvec = Vector{U}(undef, div(n * (n + 1), 2))
    MA.mutable_operate!(zero, Qvec)
    Q = SymMatrix(Qvec, n)
    for j in 1:n
        for i in 1:j
            if !iszero(Ip[j]) && !iszero(Ip[i])
                MultivariateMoments.symmetric_setindex!(
                    Q, MA.add!(Q[i, j], p.Q[Ip[i], Ip[j]]), i, j)
            end
            if !iszero(Iq[j]) && !iszero(Iq[i])
                MultivariateMoments.symmetric_setindex!(
                    Q, MA.add!(Q[i, j], q.Q[Iq[i], Iq[j]]), i, j)
            end
        end
    end
    return GramMatrix(Q, basis)
end

"""
    gram_operate(/, p::GramMatrix, α)

Computes the Gram matrix equal to the sum between `p` and `q`. On the opposite,
`p + q` gives a polynomial equal to `p + q`. The polynomial `p + q` can also be
obtained by `polynomial(gram_sum(p, q))`.
"""
function gram_operate(::typeof(/), q::GramMatrix, α)
    Q = SymMatrix(q.Q.Q / α, q.Q.n)
    return GramMatrix(Q, q.basis)
end

function Base.isapprox(p::GramMatrix, q::GramMatrix; kwargs...)
    p.basis == q.basis && isapprox(p.Q, q.Q; kwargs...)
end

struct SparseGramMatrix{T, B, U} <: AbstractGramMatrix{T, B, U}
    sub_gram_matrices::Vector{GramMatrix{T, B, U}}
end

Base.zero(::Type{SparseGramMatrix{T, B, U}}) where {T, B, U} = SparseGramMatrix(GramMatrix{T, B, U}[])
function MP.polynomial(p::SparseGramMatrix)
    return mapreduce(identity, MA.add!, p.sub_gram_matrices, init = zero(MP.polynomialtype(p)))
end
