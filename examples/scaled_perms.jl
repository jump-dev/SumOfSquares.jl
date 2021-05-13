module ScaledPerms

using GroupsCore
using Random
import AbstractAlgebra: Generic.Perm
using LinearAlgebra

abstract type MatrixGroup{T} <: Group end

abstract type MatrixGroupElem{T} <: GroupElement end

struct ScaledPermGroup{T,I} <: MatrixGroup{T}
    n::I
end

struct ScaledPerm{T,I} <: MatrixGroupElem{T}
    p::Perm{I}
    d::Diagonal{Int,Vector{T}}
end

# Group Interface

Base.one(G::ScaledPermGroup{T}) where {T} = ScaledPerm2(Perm(G.n), Diagonal(ones(T, G.n)))
GroupsCore.order(::Type{I}, G::ScaledPermGroup) = convert(I, 2^G.n * factorial(G.n))

# Base.rand
# Groups.gens

Base.eltype(::Type{ScaledPermGroup{T,I}}) where {T,I} = ScaledPerm{T,I}
# Base.IteratorSize
# Base.iterate

Base.parent(sp::ScaledPerm{T, I}) where {T, I} = ScaledPermGroup{T, I}(degree(sp.p))
GroupsCore.parent_type(::Type{ScaledPerm{T,I}}) where {T,I} = ScaledPermGroup{T,I}

Base.:(==)(sp::ScaledPerm, sq::ScaledPerm) = sp.p == sq.p && sp.d == sq.d
Base.hash(sp::ScaledPerm2, h::UInt64) = hash(sp.p, hash(sp.d, hash(ScaledPerm, u)))

function Base.deepcopy_internal(sp::ScaledPerm{T, I}, stackdict::IdDict) where {T, I}
    return ScaledPerm{T, I}(
        Base.deepcopy_internal(sp.p, stackdict),
        Base.deepcopy_internal(sp.d, stackdict)
    )
end

Base.inv(sp::ScaledPerm{T,I}) where {T,I} = ScaledPerm{T,I}(inv(sp.p), inv.(sp.d))
Base.:*(sp::ScaledPerm, sq::ScaledPerm) = ScaledPerm2(sp.p * sq.p, sp.d .* sq.d)






SymbolicWedderburn.degree(sp::ScaledPerm2) = degree(sp.p)
Base.size(sp::ScaledPerm2) = size(sp.d)

Base.@propagate_inbounds function Base.getindex(
    sp::ScaledPerm2{T},
    i::Integer,
    j::Integer,
) where {T}
    @boundscheck (1 <= i <= size(sp, 1) && 1 <= j <= size(sp, 2)) ||
                 throw(BoundsError(sp, (i, j)))
    return @inbounds p.d[i] == j ? sp.d[i, i] : zero(T)
end

Base.Matrix(m::MatrixGroupElem) = [m[i, j] for i in 1:size(m, 1) for j in 1:size(m, 2)]

LinearAlgebra.tr(sp::MatrixGroupElem) = sum(sp[i, i] for i in 1:size(sp, 1))

function SymbolicWedderburn.action_character(
    conjugacy_cls::AbstractVector{CC{<:MatrixGroupElem}},
) where CC
    vals = [tr(first(cc)) for cc in conjugacy_cls]
    return SymbolicWedderburn.Character(vals, conjugacy_cls)
end

end # of module ScaledPerms



struct ScaledPerm{T,I} <: AbstractPerm
    indices::Vector{Pair{I,T}}
end
function trace_matrix_representative(p::ScaledPerm)
    return sum(pair.second for (i, pair) in enumerate(p.indices) if pair.first == i)
end
function Base.inv(p::ScaledPerm{T,I}) where {T,I}
    indices = Vector{Pair{I,T}}(undef, length(p.indices))
    for i in eachindex(indices)
        pair = p.indices[i]
        indices[pair.first] = i => inv(pair.second)
    end
    return ScaledPerm(indices)
end

function PermutationGroups.order(p::ScaledPerm)
    cur = p
    _one = one(p)
    o = 1
    while cur != _one
        o += 1
        cur *= p
    end
    return o
end
Base.one(p::ScaledPerm{T}) where {T} =
    ScaledPerm([i => one(T) for i in eachindex(p.indices)])
Base.:(==)(p::ScaledPerm, q::ScaledPerm) = p.indices == q.indices
Base.hash(p::ScaledPerm, u::UInt64) = hash(p.indices, u)
SymbolicWedderburn.degree(p::ScaledPerm) = length(p.indices)
Base.:^(i::Integer, p::ScaledPerm) = p.indices[i]
function SymbolicWedderburn.add_inverse_permutation!(result, val, i::Int, j::Pair)
    result[i, j.first] += val / j.second
end
function SymbolicWedderburn.action_character(
    conjugacy_cls::AbstractVector{<:AbstractOrbit{<:ScaledPerm}},
)
    vals = Int[trace_matrix_representative(first(cc)) for cc in conjugacy_cls]
    return SymbolicWedderburn.Character(vals, conjugacy_cls)
end
function Base.:*(p::ScaledPerm, q::ScaledPerm)
    return ScaledPerm(map(eachindex(q.indices)) do i
        pair_q = q.indices[i]
        pair_p = p.indices[pair_q.first]
        pair_p.first => pair_p.second * pair_q.second
    end)
end
Base.:^(p::ScaledPerm, k::Integer) = Base.power_by_squaring(p, k)

function SymbolicWedderburn.permutation(
    ehom::SymbolicWedderburn.ExtensionHomomorphism{<:AbstractMonomial},
    terms::Vector{<:AbstractTerm},
)
    return ScaledPerm([ehom[monomial(term)] => coefficient(term) for term in terms])
end
