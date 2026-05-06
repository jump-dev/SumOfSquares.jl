import StarAlgebras as SA

module Monoids

import GroupsCore
import KnuthBendix as KB

import GroupsCore.gens

# /!\ Type piracy: this should go to KnuthBendix.jl
Base.convert(::Type{KB.Word{I}}, v::AbstractVector{<:Integer}) where {I} =
    KB.Word{I}(v)
Base.convert(::Type{KB.Word{I}}, w::KB.AbstractWord) where {I} = KB.Word{I}(w, false)

abstract type AbstractMonoid{I} end

function Base.iterate(M::AbstractMonoid)
    m = one(M)
    return m, (words=[word(m)], widx=1, letter=1)
end

function Base.iterate(M::AbstractMonoid, state, w = copy(state.words[state.widx]))
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

    MonoidElement{I}(
        w::AbstractVector,
        M::AbstractMonoid,
        reduced = false,
    ) where {I} = new{I,typeof(M)}(convert(KB.Word{I}, w), M, reduced)

    MonoidElement{I}(
        w::KB.Word{I},
        M::AbstractMonoid{I},
        reduced = false,
    ) where {I} = new{I,typeof(M)}(w, M, reduced)
end

function Base.deepcopy_internal(m::MonoidElement, stackdict::IdDict)
    M = parent(m)
    return M(word(m), isreduced(m))
end

# coercing to monoid
(M::AbstractMonoid{I})(
    w::AbstractVector{<:Integer},
    reduced = false,
) where {I} = MonoidElement{I}(w, M, reduced)

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

GroupsCore.gens(M::AbstractMonoid) =
    [M([i], true) for i in 1:length(KB.alphabet(M))]

# actual user-constructors for Monoid:

Base.:(/)(m::FreeMonoid{I}, rels::AbstractArray{<:MonoidElement}) where {I} =
    m / [r => one(m) for r in rels]

function Base.:(/)(
    m::FreeMonoid{I},
    rels::AbstractArray{<:Pair{<:MonoidElement,<:MonoidElement}},
    ordering = KB.LenLex
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

Base.show(io::IO, M::FreeMonoid) = print(io, "free monoid over $(KB.alphabet(M))")
Base.show(io::IO, M::Monoid) = print(
    io,
    "monoid with $(length(M.relations)) relations over $(KB.alphabet(M))",
)

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

Base.hash(m::MonoidElement, h::UInt) =
    (normalform!(m); hash(word(m), hash(parent(m), h)))

Base.show(io::IO, m::MonoidElement{I,<:FreeMonoid}) where {I} =
    print(io, KB.string_repr(word(m), KB.alphabet(parent(m))))

function Base.show(io::IO, m::MonoidElement)
    m = normalform!(m)
    return KB.print_repr(io, word(m), KB.alphabet(parent(m)), "·")
end

function wlmetric_ball(S::AbstractVector{T}; radius = 2, op = *) where {T}
    old = unique!([one(first(S)), S...])
    new = empty(old)
    sizes = [1, length(old)]
    for _ in 2:radius
        new = union!(new,
            (op(o, s) for o in @view(old[sizes[end-1]:end]) for s in S)
        )
        append!(old, new)
        resize!(new, 0)
        old = unique(old)
        push!(sizes, length(old))
    end
    return old, sizes[2:end]
end

Base.adjoint(m::MonoidElement{I}) where {I} =
    MonoidElement{I}(reverse(m.word), parent(m))

end # module

import KnuthBendix as KB

function trace_monoid(nA, nC; A=:A, C=:C)
    F, A, C = let
        SA = [Symbol(A, m) for m = 1:nA]
        SC = [Symbol(C, m) for m = 1:nC]
        F = Monoids.FreeMonoid(KB.Alphabet([SA; SC]))
        toF = Dict(zip([SA; SC], Monoids.gens(F)))

        F, [toF[e] for e in SA], [toF[e] for e in SC]
    end
    R = let
        # Reflections
        obs = [a * a => one(F) for a in [A; C]] # a² = id
        # Commutation
        com = [
            y * x => x * y for (p1, p2) in [(A, C)] for x in p1 for y in p2
        ]
        [obs; com]
    end
    M = F/R
    toM = Dict(zip([A; C], Monoids.gens(M)))
    M, [toM[e] for e in A], [toM[e] for e in C]
end

M, A, C = trace_monoid(2, 2, A=:A, C=:C)

RM = let M = M, A = A, C = C, level = 4
    A_l, sizesA = Monoids.wlmetric_ball(A, radius=level)
    C_l, sizesC = Monoids.wlmetric_ball(C, radius=level)

    # starAlg(M, 1, half = unique!([a*c for a in A_l for c in C_l]))

    @time words, sizes = Monoids.wlmetric_ball(
        unique!([a * c for a in A_l for c in C_l]);
        radius=2,
    )
    @info "Sizes of generated balls:" (A, C, combined) = (sizesA, sizesC, sizes)

    basis = SA.FixedBasis(words)
    dirac = SA.DiracMStructure(basis, *)
    table = SA.MTable(dirac, (sizes[1], sizes[1]))
    SA.StarAlgebra(M, table)
end

A = RM.(A)
C = RM.(C)
chsh = A[1] * C[1] + A[1] * C[2] + A[2] * C[1] - A[2] * C[2]

struct Full{B} <: SA.ImplicitBasis{B,B} end
Base.in(::B, ::Full{B}) where {B} = true
Base.getindex(::Full{B}, b::B) where {B} = b
import MultivariateBases as MB
MB.implicit_basis(::SA.FixedBasis{B}) where {B} = Full{B}()
MB.algebra(b::Full{B}) where {B} = SA.StarAlgebra(M, SA.DiracMStructure(b, *))
SA.comparable(::Full) = isless

f = SA.AlgebraElement(
    SA.SparseCoefficients(
        [b[k] for (k, _) in SA.nonzero_pairs(coeffs(chsh))],
        [v for (_, v) in SA.nonzero_pairs(coeffs(chsh))],
    ),
    MB.algebra(Full{eltype(SA.basis(chsh))}()),
)
n = size(b.table, 1)
gram_basis = SA.FixedBasis(SA.basis(b).elts[1:n])
one(f)
SA.coeffs(f, b)
using SumOfSquares
function SumOfSquares._term_element(a, mono::Monoids.MonoidElement)
    SA.AlgebraElement(
        SA.SparseCoefficients((mono,), (a,)),
        MB.algebra(Full{typeof(mono)}()),
    )
end

cone = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
    SA.basis(chsh),
    [gram_basis],
    [one(f)],
)
import SCS
scs = optimizer_with_attributes(
    SCS.Optimizer,
    "eps_abs" => 1e-9,
    "eps_rel" => 1e-9,
    "max_iters" => 1000,
)

import Dualization
#model = Model(Dualization.dual_optimizer(scs))
model = Model(scs)
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
MOI.Bridges.remove_bridge(backend(model).optimizer, SumOfSquares.Bridges.Constraint.ImageBridge{Float64})
@variable(model, λ)
@objective(model, Min, λ)
@constraint(model, SA.coeffs(λ * one(f) - f, b) in cone);
optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2
function _add!(f, psd, model, F, S)
    append!(psd, [
        f(MOI.get(model, MOI.ConstraintSet(), ci))
        for ci in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
    ],
    )
end
function summary(model)
    free = MOI.get(model, MOI.NumberOfVariables())
    psd = Int[]
    F = MOI.VectorOfVariables
    S = MOI.PositiveSemidefiniteConeTriangle
    _add!(MOI.side_dimension, psd, model, F, S)
    S = SCS.ScaledPSDCone
    _add!(MOI.side_dimension, psd, model, F, S)
    free -= sum(psd, init = 0) do d
        div(d * (d + 1), 2)
    end
    F = MOI.VectorAffineFunction{Float64}
    S = MOI.PositiveSemidefiniteConeTriangle
    _add!(MOI.side_dimension, psd, model, F, S)
    S = SCS.ScaledPSDCone
    _add!(MOI.side_dimension, psd, model, F, S)
    eq = Int[]
    F = MOI.VectorAffineFunction{Float64}
    S = MOI.Zeros
    _add!(MOI.dimension, eq, model, F, S)
    F = MOI.ScalarAffineFunction{Float64}
    S = MOI.EqualTo{Float64}
    _add!(MOI.dimension, eq, model, F, S)
    return free, psd, sum(eq, init = 0)
end
summary(backend(model).optimizer.model)
summary(backend(model).optimizer.model.optimizer.dual_problem.dual_model.model)
print_active_bridges(model)
