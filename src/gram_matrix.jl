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
    ::Type{S},
) where {T,B,U,S}
    return MP.polynomial_type(B, S)
end

function MP.polynomial_type(
    G::Union{AbstractGramMatrix{T,B,U},Type{<:AbstractGramMatrix{T,B,U}}},
) where {T,B,U}
    return MP.polynomial_type(G, U)
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
) where {T,B<:SA.ExplicitBasis,U}
    return GramMatrix{T,B,U,typeof(Q)}(Q, basis)
end
function GramMatrix{T,B}(
    Q::AbstractMatrix{T},
    basis::B,
) where {T,B<:SA.ExplicitBasis}
    return GramMatrix{T,B,_promote_sum(T)}(Q, basis)
end
function GramMatrix(
    # We don't use `AbstractMatrix` to avoid clash with method below with `σ`.
    Q::Union{
        MultivariateMoments.SymMatrix{T},
        MultivariateMoments.VectorizedHermitianMatrix{T},
    },
    basis::SA.ExplicitBasis,
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

SA.basis(g::GramMatrix) = g.basis
MB.implicit_basis(g::GramMatrix) = MB.implicit_basis(g.basis)

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
    basis::SA.ExplicitBasis,
    σ = 1:length(basis),
) where {T}
    return GramMatrix{T,typeof(basis)}(
        vectorized_symmetric_matrix(T, f, length(basis), σ),
        basis,
    )
end
function GramMatrix{T}(f::Function, monos::AbstractVector) where {T}
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return GramMatrix{T}(f, MB.SubBasis{MB.Monomial}(sorted_monos), σ)
end

function GramMatrix(
    Q::AbstractMatrix{T},
    basis::SA.ExplicitBasis,
    σ = 1:length(basis),
) where {T}
    return GramMatrix{T}((i, j) -> Q[σ[i], σ[j]], basis)
end
function GramMatrix(Q::AbstractMatrix, monos::AbstractVector)
    σ, sorted_monos = MP.sort_monomial_vector(monos)
    return GramMatrix(Q, MB.SubBasis{MB.Monomial}(sorted_monos), σ)
end

function MA.operate!(
    op::SA.UnsafeAddMul{typeof(*)},
    p::SA.AlgebraElement,
    g::GramMatrix,
    args::Vararg{Any,N},
) where {N}
    return MA.operate!(op, p, SA.QuadraticForm(g), args...)
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
    return show_basis_indexed_matrix(io, M)
end

struct BlockDiagonalGramMatrix{T,B,U,MT} <: AbstractGramMatrix{T,B,U}
    blocks::Vector{GramMatrix{T,B,U,MT}}
end

function MB.implicit_basis(g::BlockDiagonalGramMatrix)
    return MB.implicit_basis(first(g.blocks))
end

function MultivariateMoments.block_diagonal(blocks::Vector{<:GramMatrix})
    return BlockDiagonalGramMatrix(blocks)
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

function MA.operate!(
    op::SA.UnsafeAddMul{typeof(*)},
    p::SA.AlgebraElement,
    g::BlockDiagonalGramMatrix,
    args::Vararg{Any,N},
) where {N}
    for block in g.blocks
        MA.operate!(op, p, block, args...)
    end
    return p
end

function Base.show(io::IO, M::BlockDiagonalGramMatrix)
    print(io, "BlockDiagonalGramMatrix")
    MultivariateMoments.show_basis_indexed_blocks(io, M.blocks)
    return
end

#function Base.convert{T, PT <: AbstractPolynomial{T}}(::Type{PT}, p::GramMatrix)
#    # coefficient_type(p) may be different than T and MP.polynomial(p) may be different than PT (different module)
#    convert(PT, MP.polynomial(p))
#end

function MA.operate_to!(a::SA.AlgebraElement, ::typeof(+), g::GramMatrix)
    return MA.operate_to!(a, +, SA.QuadraticForm(g))
end

function MA.operate!(op::SA.UnsafeAdd, a::SA.AlgebraElement, g::GramMatrix)
    return MA.operate!(op, a, SA.QuadraticForm(g))
end

function MA.operate!(
    op::SA.UnsafeAdd,
    a::SA.AlgebraElement,
    g::BlockDiagonalGramMatrix,
)
    for block in g.blocks
        MA.operate!(op, a, block)
    end
    return a
end

function MA.operate_to!(
    a::SA.AlgebraElement,
    ::typeof(+),
    g::BlockDiagonalGramMatrix,
)
    MA.operate!(zero, a)
    MA.operate!(SA.UnsafeAdd(), a, g)
    MA.operate!(SA.canonical, a)
    return a
end

function MB.algebra_element(
    p::Union{GramMatrix{T,B,U},BlockDiagonalGramMatrix{T,B,U}},
) where {T,B,U}
    return MB.algebra_element(p, U)
end

function MB.algebra_element(
    g::Union{GramMatrix,BlockDiagonalGramMatrix},
    ::Type{T},
) where {T}
    a = zero(T, MB.algebra(MB.implicit_basis(g)))
    MA.operate_to!(a, +, g)
    return a
end

function MP.polynomial(
    p::Union{GramMatrix{T,B,U},BlockDiagonalGramMatrix{T,B,U}},
) where {T,B,U}
    return MP.polynomial(p, U)
end

function MP.polynomial(
    g::Union{GramMatrix,BlockDiagonalGramMatrix},
    ::Type{T},
) where {T}
    return MP.polynomial(MB.algebra_element(g, T))
end
