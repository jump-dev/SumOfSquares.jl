struct ScaledPerm{T, I} <: AbstractPerm
    indices::Vector{Pair{I, T}}
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
Base.one(p::ScaledPerm{T}) where {T} = ScaledPerm([i => one(T) for i in eachindex(p.indices)])
Base.:(==)(p::ScaledPerm, q::ScaledPerm) = p.indices == q.indices
Base.hash(p::ScaledPerm, u::UInt64) = hash(p.indices, u)
SymbolicWedderburn.degree(p::ScaledPerm) = length(p.indices)
Base.:^(i::Integer, p::ScaledPerm) = p.indices[i]
function SymbolicWedderburn.add_inverse_permutation!(result, val, i::Int, j::Pair)
    result[i, j.first] += val / j.second
end
function SymbolicWedderburn.action_character(conjugacy_cls::AbstractVector{<:AbstractOrbit{<:ScaledPerm}})
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

function SymbolicWedderburn.permutation(ehom::SymbolicWedderburn.ExtensionHomomorphism{<:AbstractMonomial}, terms::Vector{<:AbstractTerm})
    return ScaledPerm([ehom[monomial(term)] => coefficient(term) for term in terms])
end
