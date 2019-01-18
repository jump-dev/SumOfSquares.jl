export MatPolynomial

import MultivariateMoments: trimat, SymMatrix, getmat
export getmat

struct MatPolynomial{T, MT <: AbstractMonomial, MVT <: AbstractVector{MT}} <: AbstractPolynomialLike{T} # should be AbstractPolynomialLike{eltype(T)} but it doesn't work
    Q::SymMatrix{T}
    x::MVT
end
# When taking the promotion of a MatPolynomial of JuMP.Variable with a Polynomial JuMP.Variable, it should be a Polynomial of AffExpr
MP.coefficienttype(::Type{<:MatPolynomial{T}}) where {T} = Base.promote_op(+, T, T)
MP.polynomialtype(::Type{MatPolynomial{T, MT, MVT}}) where {T, MT, MVT} = polynomialtype(MT, coefficienttype(MatPolynomial{T, MT, MVT}))
MP.polynomialtype(::Type{MatPolynomial{T, MT, MVT}}, ::Type{S}) where {S, T, MT, MVT} = polynomialtype(MT, S)
MP.variables(p::MatPolynomial) = MP.variables(p.x)
MP.nvariables(p::MatPolynomial) = MP.nvariables(p.x)

Base.zero(::Type{MatPolynomial{T, MT, MVT}}) where {T, MT, MVT} = MatPolynomial{T, MT, monovectype(MT)}(SymMatrix{T}(T[], 0), emptymonovec(MT))
Base.iszero(p::MatPolynomial) = iszero(polynomial(p))

Base.getindex(p::MatPolynomial, I...) = getindex(p.Q, I...)
Base.copy(p::MatPolynomial) = MatPolynomial(copy(p.Q), copy(p.x))

MultivariateMoments.getmat(p::MatPolynomial{T}) where {T} = p.Q

function _matplus(p::MatPolynomial, q::MatPolynomial)
    @assert p.x == q.x
    MatPolynomial(p.Q+q.Q, p.x)
end

function MatPolynomial{T}(f::Function, x::AbstractVector{MT}, σ) where {T, MT<:AbstractMonomial}
    MatPolynomial{T, MT, monovectype(x)}(trimat(T, f, length(x), σ), x)
end
MatPolynomial{T}(f::Function, x::AbstractVector, σ) where T = MatPolynomial{T}(f, monomial.(x), σ)
function MatPolynomial{T}(f::Function, x::AbstractVector) where T
    σ, X = sortmonovec(x)
    MatPolynomial{T}(f, X, σ)
end
MatPolynomial(f::Function, x) = MatPolynomial{Base.promote_op(f, Int, Int)}(f, x)

function MatPolynomial(Q::AbstractMatrix{T}, x, σ) where T
    MatPolynomial{T}((i,j) -> Q[σ[i], σ[j]], x)
end
function MatPolynomial(Q::AbstractMatrix, x)
    σ, X = sortmonovec(x)
    MatPolynomial(Q, X, σ)
end

#function Base.convert{T, PT <: AbstractPolynomial{T}}(::Type{PT}, p::MatPolynomial)
#    # coefficienttype(p) may be different than T and polynomial(p) may be different than PT (different module)
#    convert(PT, polynomial(p))
#end
function MP.polynomial(p::MatPolynomial)
    MP.polynomial(getmat(p), p.x)
end
function MP.polynomial(p::MatPolynomial, ::Type{S}) where {S}
    MP.polynomial(getmat(p), p.x, S)
end

const APL = AbstractPolynomialLike

Base.:(+)(x::APL, y::MatPolynomial) = x + polynomial(y)
Base.:(+)(x::MatPolynomial, y::APL) = polynomial(x) + y
Base.:(+)(x::MatPolynomial, y::MatPolynomial) = polynomial(x) + polynomial(y)
Base.:(-)(x::APL, y::MatPolynomial) = x - polynomial(y)
Base.:(-)(x::MatPolynomial, y::APL) = polynomial(x) - y
Base.:(-)(x::MatPolynomial, y::MatPolynomial) = polynomial(x) - polynomial(y)
Base.:(==)(p::APL, q::MatPolynomial) = p == polynomial(q)
Base.:(==)(p::MatPolynomial, q::MatPolynomial) = iszero(p - q)
function Base.isapprox(p::MatPolynomial, q::MatPolynomial; kwargs...)
    p.x == q.x && isapprox(p.Q, q.Q; kwargs...)
end
