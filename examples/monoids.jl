# MIT License
# 
# Copyright (c) 2022 Marek Kaluba <kalmar@mailbox.org>
#
# This `Monoids` module and `trace_monoid` function are adapted
# from a private package of Marek Kaluba.

module Monoids

import GroupsCore
import KnuthBendix as KB

# /!\ Type piracy: this should go to KnuthBendix.jl
function Base.convert(
    ::Type{KB.Word{I}},
    v::AbstractVector{<:Integer},
) where {I}
    return KB.Word{I}(v)
end
function Base.convert(::Type{KB.Word{I}}, w::KB.AbstractWord) where {I}
    return KB.Word{I}(w, false)
end

abstract type AbstractMonoid{I} end

function Base.iterate(M::AbstractMonoid)
    m = one(M)
    return m, (words=[word(m)], widx=1, letter=1)
end

function Base.iterate(
    M::AbstractMonoid,
    state,
    w=copy(state.words[state.widx]),
)
    next_letter = state.letter ≥ length(KB.alphabet(M)) ? 1 : state.letter + 1
    next_widx = next_letter == 1 ? state.widx + 1 : state.widx
    next_state = (words=state.words, letter=next_letter, widx=next_widx)

    push!(w, state.letter)
    m = M(w) # w gets reduced here
    normalform!(m)
    if word(m) in state.words
        w = word(m)
        KB.Words.store!(w, state.words[state.widx])
        Base.iterate(M, next_state, w) # we pass w for reuse
    else
        push!(next_state.words, word(m))
        return m, next_state # we'll need to allocate new w
    end
end

Base.IteratorSize(::Type{<:AbstractMonoid}) = Base.SizeUnknown()

Base.eltype(::Type{M}) where {I,M<:AbstractMonoid{I}} = MonoidElement{I,M}

# types and constructors

struct FreeMonoid{I,A} <: AbstractMonoid{I}
    alphabet::A

    function FreeMonoid{I}(a::KB.Alphabet) where {I<:Integer}
        length(a) < typemax(UInt16) ||
            throw("Too many letters in alphabet. Try with $(widen(I))")
        return new{I,typeof(a)}(a)
    end
end

FreeMonoid(a::KB.Alphabet) = FreeMonoid{UInt16}(a)
FreeMonoid(n::Integer) = FreeMonoid(KB.Alphabet([Symbol('a', i) for i in 1:n]))

const Relation{I} = Pair{KB.Word{I},KB.Word{I}}

struct Monoid{I,A,R} <: AbstractMonoid{I}
    alphabet::A
    relations::Vector{Relation{I}} # in case you need them later
    rws::R
end

mutable struct MonoidElement{I,M}
    word::KB.Word{I}
    parent::M
    reduced::Bool

    function MonoidElement{I}(
        w::AbstractVector,
        M::AbstractMonoid,
        reduced=false,
    ) where {I}
        return new{I,typeof(M)}(convert(KB.Word{I}, w), M, reduced)
    end

    function MonoidElement{I}(
        w::KB.Word{I},
        M::AbstractMonoid{I},
        reduced=false,
    ) where {I}
        return new{I,typeof(M)}(w, M, reduced)
    end
end

function Base.deepcopy_internal(m::MonoidElement, stackdict::IdDict)
    M = parent(m)
    return M(word(m), isreduced(m))
end

# coercing to monoid
function (M::AbstractMonoid{I})(
    w::AbstractVector{<:Integer},
    reduced=false,
) where {I}
    return MonoidElement{I}(w, M, reduced)
end

# Accessors and basic manipulation
KB.alphabet(M::FreeMonoid) = M.alphabet
KB.alphabet(M::Monoid) = M.alphabet
rewriting(M::Monoid) = M.rws

Base.parent(m::MonoidElement) = m.parent
word(m::MonoidElement) = m.word
isreduced(m::MonoidElement) = m.reduced

Base.one(M::AbstractMonoid{I}) where {I} = M(one(KB.Word{I}), true)
Base.one(m::MonoidElement) = one(parent(m))
# isone is equivalent to the word problem for general inputs, see also ==
# isone(m::MonoidElement) = one(parent(m)) == m # should be literal empty word comparison, but this is already implemented in Base.:(==) below

function GroupsCore.gens(M::AbstractMonoid)
    return [M([i], true) for i in 1:length(KB.alphabet(M))]
end

# actual user-constructors for Monoid:

function Base.:(/)(
    m::FreeMonoid{I},
    rels::AbstractArray{<:MonoidElement},
) where {I}
    return m / [r => one(m) for r in rels]
end

function Base.:(/)(
    m::FreeMonoid{I},
    rels::AbstractArray{<:Pair{<:MonoidElement,<:MonoidElement}},
    ordering=KB.LenLex,
) where {I}
    A = m.alphabet
    new_rels = Relation{I}[word(first(r)) => word(last(r)) for r in rels]

    rws = KB.RewritingSystem(
        Tuple{KB.Word{I},KB.Word{I}}[(first(r), last(r)) for r in new_rels],
        ordering(KB.alphabet(m)),
    )

    rws = KB.knuthbendix(KB.Settings(), rws)

    return Monoid(A, new_rels, rws)
end

function Base.show(io::IO, M::FreeMonoid)
    return print(io, "free monoid over $(KB.alphabet(M))")
end
function Base.show(io::IO, M::Monoid)
    return print(
        io,
        "monoid with $(length(M.relations)) relations over $(KB.alphabet(M))",
    )
end

function Base.:(*)(m1::MonoidElement, m2::MonoidElement)
    parent(m1) === parent(m2) ||
        throw("cannot multiply elements from different monoids")
    return parent(m1)(word(m1) * word(m2))
end

function Base.in(m1::MonoidElement, m2::MonoidElement)
    parent(m1) === parent(m2) ||
        throw("cannot compare elements from different monoids")
    return issubset(word(m1), word(m2))
end

Base.:(^)(m::MonoidElement, n::Integer) = Base.power_by_squaring(m, n)

"""
    normalform!(m::MonoidElement[, tmp::AbstractWord])
Reduce `m` to its normalform, as defined by the rewriting of `parent(m)`.
"""
function normalform!(m::MonoidElement)
    return isreduced(m) ? m : normalform!(m, similar(m.word))
end

function normalform!(m::MonoidElement, tmp::KB.AbstractWord)
    if !isreduced(m)
        resize!(tmp, length(m.word))
        copyto!(tmp, m.word)
        resize!(m.word, 0)
        KB.rewrite!(m.word, tmp, rewriting(parent(m)))
        m.reduced = true
    end
    @assert isreduced(m)
    return m
end

function Base.:(==)(m1::MonoidElement, m2::MonoidElement)
    parent(m1) === parent(m2) || return false
    word(m1) == word(m2) && return true
    normalform!(m1)
    normalform!(m2)
    return word(m1) == word(m2)
end

function Base.isless(m1::MonoidElement, m2::MonoidElement)
    parent(m1) === parent(m2) ||
        throw("cannot compare elements from different monoids")
    normalform!(m1)
    normalform!(m2)
    return isless(word(m1), word(m2))
end

function Base.hash(m::MonoidElement, h::UInt)
    return (normalform!(m); hash(word(m), hash(parent(m), h)))
end

function Base.show(io::IO, m::MonoidElement{I,<:FreeMonoid}) where {I}
    return print(io, KB.string_repr(word(m), KB.alphabet(parent(m))))
end

function Base.show(io::IO, m::MonoidElement)
    m = normalform!(m)
    return KB.print_repr(io, word(m), KB.alphabet(parent(m)), "·")
end

function wlmetric_ball(S::AbstractVector{T}; radius=2, op=*) where {T}
    old = unique!([one(first(S)), S...])
    new = empty(old)
    sizes = [1, length(old)]
    for _ in 2:radius
        new = union!(
            new,
            (op(o, s) for o in @view(old[sizes[end-1]:end]) for s in S),
        )
        append!(old, new)
        resize!(new, 0)
        old = unique(old)
        push!(sizes, length(old))
    end
    return old, sizes[2:end]
end

function Base.adjoint(m::MonoidElement{I}) where {I}
    return MonoidElement{I}(reverse(m.word), parent(m))
end

end # module

function trace_monoid(nA, nC; A=:A, C=:C)
    F, A, C = let
        SA = [Symbol(A, m) for m in 1:nA]
        SC = [Symbol(C, m) for m in 1:nC]
        F = Monoids.FreeMonoid(KB.Alphabet([SA; SC]))
        toF = Dict(zip([SA; SC], GroupsCore.gens(F)))

        F, [toF[e] for e in SA], [toF[e] for e in SC]
    end
    R = let
        # Reflections
        obs = [a * a => one(F) for a in [A; C]] # a² = id
        # Commutation
        com = [y * x => x * y for (p1, p2) in [(A, C)] for x in p1 for y in p2]
        [obs; com]
    end
    M = F / R
    toM = Dict(zip([A; C], GroupsCore.gens(M)))
    return M, [toM[e] for e in A], [toM[e] for e in C]
end
