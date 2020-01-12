export GramMatrix, SparseGramMatrix

import MultivariateMoments: trimat, SymMatrix, getmat
export gram_operate, getmat

abstract type AbstractGramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}, U} <: MP.AbstractPolynomialLike{U} end

function MP.monomialtype(::Union{AbstractGramMatrix{T, MT},
                                 Type{<:AbstractGramMatrix{T, MT}}}) where {T, MT}
    return MT
end
MP.polynomialtype(::Union{AbstractGramMatrix{T, MT, MVT, U}, Type{<:AbstractGramMatrix{T, MT, MVT, U}}}) where {T, MT, MVT, U} = MP.polynomialtype(MT, U)

Base.:(+)(x::MP.APL, y::AbstractGramMatrix) = x + MP.polynomial(y)
Base.:(+)(x::AbstractGramMatrix, y::MP.APL) = MP.polynomial(x) + y
Base.:(+)(x::AbstractGramMatrix, y::AbstractGramMatrix) = MP.polynomial(x) + MP.polynomial(y)
Base.:(-)(x::MP.APL, y::AbstractGramMatrix) = x - MP.polynomial(y)
Base.:(-)(x::AbstractGramMatrix, y::MP.APL) = MP.polynomial(x) - y
Base.:(-)(x::AbstractGramMatrix, y::AbstractGramMatrix) = MP.polynomial(x) - MP.polynomial(y)
Base.:(==)(p::MP.APL, q::AbstractGramMatrix) = p == MP.polynomial(q)
Base.:(==)(p::AbstractGramMatrix, q::AbstractGramMatrix) = iszero(p - q)

"""
    struct GramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} <: MP.APL{T}
        Q::SymMatrix{T}
        x::MVT
    end

Gram matrix ``x^\\top Q x`` where `Q` is a symmetric matrix indexed by the
vector of monomials `x`.
"""
struct GramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}, U} <: AbstractGramMatrix{T, MT, MVT, U}
    Q::SymMatrix{T}
    x::MVT
end
GramMatrix{T, MT, MVT}(Q::SymMatrix{T}, x::MVT) where {T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} = GramMatrix{T, MT, MVT, _promote_sum(T)}(Q, x)
GramMatrix(Q::SymMatrix{T}, x::MVT) where {T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} = GramMatrix{T, MT, MVT}(Q, x)
# When taking the promotion of a GramMatrix of JuMP.Variable with a Polynomial JuMP.Variable, it should be a Polynomial of AffExpr
function MP.constantmonomial(p::GramMatrix)
    if isempty(p.x)
        return MP.constantmonomial(MP.monomialtype(p))
    else
        return MP.constantmonomial(p.x[1])
    end
end
MP.variables(p::GramMatrix) = MP.variables(p.x)
MP.nvariables(p::GramMatrix) = MP.nvariables(p.x)

Base.zero(::Type{GramMatrix{T, MT, MVT, U}}) where {T, MT, MVT, U} = GramMatrix{T, MT, MVT, U}(SymMatrix{T}(T[], 0), MP.emptymonovec(MT))
Base.iszero(p::GramMatrix) = iszero(MP.polynomial(p))

Base.getindex(p::GramMatrix, I...) = getindex(p.Q, I...)
Base.copy(p::GramMatrix) = GramMatrix(copy(p.Q), copy(p.x))

MultivariateMoments.getmat(p::GramMatrix{T}) where {T} = p.Q

function GramMatrix{T}(f::Function, x::AbstractVector{MT}, σ) where {T, MT<:MP.AbstractMonomial}
    GramMatrix{T, MT, MP.monovectype(x)}(trimat(T, f, length(x), σ), x)
end
GramMatrix{T}(f::Function, x::AbstractVector, σ) where T = GramMatrix{T}(f, MP.monomial.(x), σ)
function GramMatrix{T}(f::Function, x::AbstractVector) where T
    σ, X = MP.sortmonovec(x)
    GramMatrix{T}(f, X, σ)
end
GramMatrix(f::Function, x) = GramMatrix{Base.promote_op(f, Int, Int)}(f, x)

function GramMatrix(Q::AbstractMatrix{T}, x, σ) where T
    GramMatrix{T}((i,j) -> Q[σ[i], σ[j]], x)
end
function GramMatrix(Q::AbstractMatrix, x)
    σ, X = MP.sortmonovec(x)
    GramMatrix(Q, X, σ)
end

#function Base.convert{T, PT <: AbstractPolynomial{T}}(::Type{PT}, p::GramMatrix)
#    # coefficienttype(p) may be different than T and MP.polynomial(p) may be different than PT (different module)
#    convert(PT, MP.polynomial(p))
#end
function MP.polynomial(p::GramMatrix{T, MT, MVT, U}) where {T, MT, MVT, U}
    return MP.polynomial(p, U)
end
function MP.polynomial(p::GramMatrix, ::Type{S}) where {S}
    return MP.polynomial(getmat(p), p.x, S)
end

# The `i`th index of output is the index of occurence of `x[i]` in `y`,
# or `0` if it does not occur.
function multi_findsorted(x, y)
    I = zeros(Int, length(x))
    j = 1
    for i in eachindex(x)
        while j ≤ length(y) && x[i] < y[j]
            j += 1
        end
        if j ≤ length(y) && x[i] == y[j]
            I[i] = j
        end
    end
    return I
end

"""
    gram_operate(::typeof(+), p::GramMatrix{S}, q::GramMatrix{T})

Computes the Gram matrix equal to the sum between `p` and `q`. On the opposite,
`p + q` gives a polynomial equal to `p + q`. The polynomial `p + q` can also be
obtained by `polynomial(gram_sum(p, q))`.
"""
function gram_operate(::typeof(+), p::GramMatrix{S}, q::GramMatrix{T}) where {S, T}
    monos = MultivariatePolynomials.mergemonovec([p.x, q.x])
    U = typeof(zero(S) + zero(T))
    n = length(monos)
    Q = SymMatrix(zeros(U, div(n * (n + 1), 2)), n)
    Ip = multi_findsorted(monos, p.x)
    Iq = multi_findsorted(monos, q.x)
    for j in 1:length(monos)
        for i in 1:j
            if !iszero(Ip[j]) && !iszero(Ip[i])
                MultivariateMoments.symmetric_setindex!(
                    Q, Q[i, j] + p.Q[Ip[i], Ip[j]], i, j)
            end
            if !iszero(Iq[j]) && !iszero(Iq[i])
                MultivariateMoments.symmetric_setindex!(
                    Q, Q[i, j] + q.Q[Iq[i], Iq[j]], i, j)
            end
        end
    end
    return GramMatrix(Q, monos)
end

"""
    gram_operate(/, p::GramMatrix{S}, q::GramMatrix{T})

Computes the Gram matrix equal to the sum between `p` and `q`. On the opposite,
`p + q` gives a polynomial equal to `p + q`. The polynomial `p + q` can also be
obtained by `polynomial(gram_sum(p, q))`.
"""
function gram_operate(::typeof(/), q::GramMatrix, α)
    Q = SymMatrix(q.Q.Q / α, q.Q.n)
    return GramMatrix(Q, q.x)
end

function Base.isapprox(p::GramMatrix, q::GramMatrix; kwargs...)
    p.x == q.x && isapprox(p.Q, q.Q; kwargs...)
end

struct SparseGramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}, U} <: AbstractGramMatrix{T, MT, MVT, U}
    sub_gram_matrices::Vector{GramMatrix{T, MT, MVT, U}}
end

Base.zero(::Type{SparseGramMatrix{T, MT, MVT, U}}) where {T, MT, MVT, U} = SparseGramMatrix(GramMatrix{T, MT, MVT, U}[])
function MP.polynomial(p::SparseGramMatrix)
    return mapreduce(identity, MA.add!, p.sub_gram_matrices, init = zero(MP.polynomialtype(p)))
end
