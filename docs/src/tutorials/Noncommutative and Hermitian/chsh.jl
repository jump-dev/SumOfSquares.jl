# # CHSH

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Noncommutative and Hermitian/chsh.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Noncommutative and Hermitian/chsh.ipynb)
# **Contributed by**: Marek Kaluba and Benoît Legat
# **Adapted from**: [Talk](https://jump.dev/assets/jump-dev-workshops/2024/legat.html) at [JuMP-dev 2024](https://jump.dev/meetings/jumpdev2024/)
#
# The goal of this tutorial is to show how to use SumOfSquares.jl with a custom algebra that is **not** defined
# with [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl/) or
# [TypedPolynomials](https://github.com/JuliaAlgebra/TypedPolynomials.jl/).
# Even though some part of SumOfSquares.jl only works for monomials defined by these packages,
# there is an effort to abstract away as much as possible on top of
# [StarAlgebras](https://github.com/JuliaAlgebra/StarAlgebras.jl/).
#
# This illustrate this, in this tutorial, we are doing to define a custom monoid structure implementing the
# rewriting rules of the [CHSH inequality](https://en.wikipedia.org/wiki/CHSH_inequality) using
# [KnuthBendix](https://github.com/kalmarek/KnuthBendix.jl).
# So we start with this `Monoids` module and `trace_monoid` function adapted
# from internal notes of Marek Kaluba.

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

import StarAlgebras as SA
import KnuthBendix as KB
import GroupsCore

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

# The rewriting rule are as follows:

monoid, A, C = trace_monoid(2, 2, A=:A, C=:C)
monoid.rws

# We now define a `StarAlgebra` from [StarAlgebras](https://github.com/JuliaAlgebra/StarAlgebras.jl/)

RM = let monoid = monoid, A = A, C = C, level = 4
    A_l, sizesA = Monoids.wlmetric_ball(A, radius=level)
    C_l, sizesC = Monoids.wlmetric_ball(C, radius=level)

    # starAlg(M, 1, half = unique!([a*c for a in A_l for c in C_l]))

    @time words, sizes = Monoids.wlmetric_ball(
        unique!([a * c for a in A_l for c in C_l]);
        radius=2,
    )
    @info "Sizes of generated balls:" (A, C, combined) =
        (sizesA, sizesC, sizes)

    basis = SA.FixedBasis(words)
    dirac = SA.DiracMStructure(basis, *)
    table = SA.MTable(dirac, (sizes[1], sizes[1]))
    SA.StarAlgebra(monoid, table)
end

# We can convert a monoid element:

A[1], typeof(A[1])

# to an element of the algebra `RM` (so essentially `1 ⋅ A`) as follows:

RM(A[1])

# Then, we can do arithmetic on these algebra element to form the CHSH inequality:

chsh = let A = RM.(A), C = RM.(C)
    A[1] * C[1] + A[1] * C[2] + A[2] * C[1] - A[2] * C[2]
end

# SumOfSquares needs the ambient implicit basis containing all the monoid elements
# so we define it as follows:

struct Full{B} <: SA.ImplicitBasis{B,B} end
Base.in(::B, ::Full{B}) where {B} = true
Base.getindex(::Full{B}, b::B) where {B} = b
import MultivariateBases as MB
MB.implicit_basis(::SA.FixedBasis{B}) where {B} = Full{B}()
function MB.algebra(b::Full{B}) where {B}
    return SA.StarAlgebra(monoid, SA.DiracMStructure(b, *))
end
SA.comparable(::Full) = isless

# The `chsh` polynomial can be rewritten in this basis as follows:

f = SA.AlgebraElement(
    SA.SparseCoefficients(
        [SA.basis(chsh)[k] for (k, _) in SA.nonzero_pairs(SA.coeffs(chsh))],
        [v for (_, v) in SA.nonzero_pairs(SA.coeffs(chsh))],
    ),
    MB.algebra(Full{eltype(SA.basis(chsh))}()),
)

# We pick the SCS solver:

using SumOfSquares
import SCS
scs = optimizer_with_attributes(
    SCS.Optimizer,
    "eps_abs" => 1e-9,
    "eps_rel" => 1e-9,
    "max_iters" => 1000,
)

# We currently need to manually add all SumOfSquares bridges as follows:

model = Model(scs)
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
@variable(model, λ)
@objective(model, Min, λ)

# Adding the SumOfSquares constraint is currently not as user-friendly than
# with monomials:

n = size(SA.mstructure(RM).table, 1)
gram_basis = SA.FixedBasis(SA.basis(chsh).elts[1:n])
cone = SumOfSquares.WeightedSOSCone{MOI.PositiveSemidefiniteConeTriangle}(
    SA.basis(chsh),
    [gram_basis],
    [one(f)],
)
@constraint(model, SA.coeffs(λ * one(f) - f, SA.basis(chsh)) in cone);

optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2

# Let's look at the size of the generated SDP:

function _add!(f, psd, model, F, S)
    return append!(
        psd,
        [
            f(MOI.get(model, MOI.ConstraintSet(), ci)) for
            ci in MOI.get(model, MOI.ListOfConstraintIndices{F,S}())
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
    free -= sum(psd, init=0) do d
        return div(d * (d + 1), 2)
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
    println(
        "$free free variables, $(sum(eq, init = 0)) equality constraints, PSD block sizes: $psd",
    )
    return
end
summary(backend(model).optimizer.model)

# We can see that the `FreeBridge` was used because SCS supports free variables and
# affine-in-PSD constraints.

print_active_bridges(model)

# As a workaround for [this MOI issue](https://github.com/jump-dev/MathOptInterface.jl/pull/3001),
# we need to remove the Image bridge

import Dualization
model = Model(Dualization.dual_optimizer(scs))
SumOfSquares.Bridges.add_all_bridges(backend(model).optimizer, Float64)
MOI.Bridges.remove_bridge(
    backend(model).optimizer,
    SumOfSquares.Bridges.Constraint.ImageBridge{Float64},
)
@variable(model, λ)
@objective(model, Min, λ)
@constraint(model, SA.coeffs(λ * one(f) - f, SA.basis(chsh)) in cone);
optimize!(model)
solution_summary(model)
objective_value(model) ≈ 2√2

# The model is much smaller this time, with 289 equality constraints and 1 free variable.

summary(backend(model).optimizer.model)

# After dualization that is 289 free variables and 1 equality constraints.
# This is a big improvement compared to the 3051 free variables, 18 equality constraints without dualization.

summary(backend(model).optimizer.model.optimizer.dual_problem.dual_model)

# We can see that the bridge used is the `KernelBridge` this time.

print_active_bridges(model)

# We didn't test using the `ImageBridge` with dualization or using the `KernelBridge` without dualization.
# However, these will always produce larger models, which is why the choice o bridge can be need
# automatically once you choose whether you want to dualize the problem or not.
