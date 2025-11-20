using PermutationGroups
d = perm"(1, 2, 3, 4)"
c = perm"(1, 3)"
G = PermGroup([c, d])

import GroupsCore

struct DihedralGroup <: GroupsCore.Group
    n::Int
end

struct DihedralElement <: GroupsCore.GroupElement
    n::Int
    reflection::Bool
    id::Int
end

Base.one(G::DihedralGroup) = DihedralElement(G.n, false, 0)

Base.eltype(::Type{DihedralGroup}) = DihedralElement
Base.IteratorSize(::Type{DihedralGroup}) = Base.HasLength()

function Base.iterate(G::DihedralGroup, prev::DihedralElement=DihedralElement(G.n, false, -1))
    if prev.id + 1 >= G.n
        if prev.reflection
            return nothing
        else
            next = DihedralElement(G.n, true, 0)
        end
    else
        next = DihedralElement(G.n, prev.reflection, prev.id + 1)
    end
    return next, next
end

GroupsCore.order(::Type{T}, G::DihedralGroup) where {T} = convert(T, 2G.n)
GroupsCore.gens(G::DihedralGroup) = [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]

Base.parent(g::DihedralElement) = DihedralGroup(g.n)
function Base.:(==)(g::DihedralElement, h::DihedralElement)
    return g.n == h.n && g.reflection == h.reflection && g.id == h.id
end

function Base.inv(el::DihedralElement)
    if el.reflection || iszero(el.id)
        return el
    else
        return DihedralElement(el.n, false, el.n - el.id)
    end
end
function Base.:*(a::DihedralElement, b::DihedralElement)
    a.n == b.n || error("Cannot multiply elements from different Dihedral groups")
    id = mod(a.reflection ? a.id - b.id : a.id + b.id, a.n)
    return DihedralElement(a.n, a.reflection != b.reflection, id)
end

Base.copy(a::DihedralElement) = DihedralElement(a.n, a.reflection, a.id)

function GroupsCore.order(el::DihedralElement)
    if el.reflection
        return 2
    else
        if iszero(el.id)
            return 1
        else
            return div(el.n, gcd(el.n, el.id))
        end
    end
end

using SumOfSquares
using DynamicPolynomials
@polyvar x y
struct DihedralAction <: Symmetry.OnMonomials end
import SymbolicWedderburn
SymbolicWedderburn.coeff_type(::DihedralAction) = Float64
function SymbolicWedderburn.action(::DihedralAction, el::DihedralElement, mono::AbstractMonomial)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return mono([x, y] => [sign_x * var_x, sign_y * var_y])
end

poly = x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1

G = DihedralGroup(4)
for g in G
    @show SymbolicWedderburn.action(DihedralAction(), g, poly)
end

import CSDP
function solve(G)
    solver = CSDP.Optimizer
    model = Model(solver)
    @variable(model, t)
    @objective(model, Max, t)
    pattern = Symmetry.Pattern(G, DihedralAction())
    con_ref = @constraint(model, poly - t in SOSCone(), symmetry = pattern)
    optimize!(model)
    @show value(t)


    for g in gram_matrix(con_ref).blocks
        println(g.basis.polynomials)
    end
end
solve(G)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
