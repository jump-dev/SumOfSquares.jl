export GramMatrix

import MultivariateMoments: trimat, SymMatrix, getmat
export getmat

"""
    struct GramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} <: MP.APL{T}
        Q::SymMatrix{T}
        x::MVT
    end

Gram matrix ``x^\\top Q x`` where `Q` is a symmetric matrix indexed by the
vector of monomials `x`.
"""
struct GramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} <: MP.APL{T} # should be MP.APL{eltype(T)} but it doesn't work
    Q::SymMatrix{T}
    x::MVT
end
# When taking the promotion of a GramMatrix of JuMP.Variable with a Polynomial JuMP.Variable, it should be a Polynomial of AffExpr
MP.coefficienttype(::Type{<:GramMatrix{T}}) where {T} = Base.promote_op(+, T, T)
function MP.constantmonomial(p::GramMatrix)
    if isempty(p.x)
        return MP.constantmonomial(MP.monomialtype(p))
    else
        return MP.constantmonomial(p.x[1])
    end
end
function MP.monomialtype(::Union{GramMatrix{T, MT},
                                 Type{<:GramMatrix{T, MT}}}) where {T, MT}
    return MT
end
MP.polynomialtype(::Type{GramMatrix{T, MT, MVT}}) where {T, MT, MVT} = MP.polynomialtype(MT, MP.coefficienttype(GramMatrix{T, MT, MVT}))
MP.polynomialtype(::Type{GramMatrix{T, MT, MVT}}, ::Type{S}) where {S, T, MT, MVT} = MP.polynomialtype(MT, S)
MP.variables(p::GramMatrix) = MP.variables(p.x)
MP.nvariables(p::GramMatrix) = MP.nvariables(p.x)

Base.zero(::Type{GramMatrix{T, MT, MVT}}) where {T, MT, MVT} = GramMatrix{T, MT, MP.monovectype(MT)}(SymMatrix{T}(T[], 0), MP.emptymonovec(MT))
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
function MP.polynomial(p::GramMatrix)
    MP.polynomial(getmat(p), p.x)
end
function MP.polynomial(p::GramMatrix, ::Type{S}) where {S}
    MP.polynomial(getmat(p), p.x, S)
end

Base.:(+)(x::MP.APL, y::GramMatrix) = x + MP.polynomial(y)
Base.:(+)(x::GramMatrix, y::MP.APL) = MP.polynomial(x) + y
Base.:(+)(x::GramMatrix, y::GramMatrix) = MP.polynomial(x) + MP.polynomial(y)
Base.:(-)(x::MP.APL, y::GramMatrix) = x - MP.polynomial(y)
Base.:(-)(x::GramMatrix, y::MP.APL) = MP.polynomial(x) - y
Base.:(-)(x::GramMatrix, y::GramMatrix) = MP.polynomial(x) - MP.polynomial(y)
Base.:(==)(p::MP.APL, q::GramMatrix) = p == MP.polynomial(q)
Base.:(==)(p::GramMatrix, q::GramMatrix) = iszero(p - q)
function Base.isapprox(p::GramMatrix, q::GramMatrix; kwargs...)
    p.x == q.x && isapprox(p.Q, q.Q; kwargs...)
end
