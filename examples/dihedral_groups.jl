module DihedralGroups

using GroupsCore
using Random

export DihedralGroup, DihedralElement

struct DihedralGroup <: Group
    n::Int
end

struct DihedralElement <: GroupElement
    reflection::Bool
    id::Int
    parent::DihedralGroup
    DihedralElement(refl::Bool, id::Integer, D::DihedralGroup) = new(refl, mod(id, D.n), D)
end

isreflection(el::DihedralElement) = el.reflection

# Group Interface:
Base.one(D::DihedralGroup) = DihedralElement(false, 0, D)
GroupsCore.order(::Type{T}, D::DihedralGroup) where {T} = convert(T, 2D.n)
GroupsCore.gens(G::DihedralGroup) =
    [DihedralElement(false, 1, G), DihedralElement(true, 0, G)]

Base.rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{DihedralGroup}) =
    (D = rs[]; DihedralElement(rand(rng, Bool), rand(rng, 0:D.n-1), D))

Base.eltype(::Type{DihedralGroup}) = DihedralElement
Base.IteratorSize(::Type{DihedralGroup}) = Base.HasLength()

function Base.iterate(D::DihedralGroup, state = 0)
    state >= order(D) && return nothing
    elt =
        iseven(state) ? DihedralElement(false, state ÷ 2, D) :
        DihedralElement(true, (state - 1) ÷ 2, D)
    return elt, state + 1
end

# GroupElement interface
Base.parent(el::DihedralElement) = el.parent
GroupsCore.parent_type(::Type{DihedralElement}) = DihedralGroup

Base.:(==)(g::DihedralElement, h::DihedralElement) =
    (parent(g) === parent(h) && isreflection(g) == isreflection(h) && g.id == h.id)

function Base.inv(el::DihedralElement)
    (isone(el) || isreflection(el)) && return el
    D = parent(el)
    return DihedralElement(false, D.n - el.id, D)
end

function Base.:*(a::DihedralElement, b::DihedralElement)
    parent(a) === parent(b) ||
        error("Cannot multiply elements from different Dihedral groups")
    D = parent(a)
    id = isreflection(a) ? a.id - b.id : a.id + b.id
    return DihedralElement(isreflection(a) != isreflection(b), id, D)
end

# Not necessary for Group Interface:
Base.show(io::IO, D::DihedralGroup) = print(io, "dihedral group of $(D.n)-gon")
Base.show(io::IO, g::DihedralElement) = print(io, (isreflection(g) ? "R" : "C"), g.id)


#=
# testing
using Test
include(joinpath(dirname(pathof(GroupsCore)), "..", "test", "conformance_test.jl"))

D = DihedralGroup(6)
test_Group_interface(D)
test_GroupElement_interface(rand(D, 2)...)
=#

end # of module DihedralGroups
