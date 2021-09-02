# # Dihedral symmetry of the Robinson form

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry/dihedral_symmetry_of_the_robinson_form.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry/dihedral_symmetry_of_the_robinson_form.ipynb)
# **Adapted from**: Example 5.4 of [GP04]
#
# [GP04] Gatermann, Karin and Parrilo, Pablo A.
# *Symmetry groups, semidefinite programs, and sums of squares*.
# Journal of Pure and Applied Algebra 192.1-3 (2004): 95-128.


using Test #src

# We start by defining the Dihedral group of order 8.
# This group is isomorphic to the following permutation group:

using PermutationGroups
d = perm"(1, 2, 3, 4)"
c = perm"(1, 3)"
G = PermGroup([c, d])

# We could rely on this isomorphism to define this group.
# However, in order to illustrate how to do symmetry reduction with a custom group,
# we show in this example what should be implemented to define a new group.

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

Base.eltype(::DihedralGroup) = DihedralElement
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
Base.IteratorSize(::Type{DihedralGroup}) = Base.HasLength()

GroupsCore.order(::Type{T}, G::DihedralGroup) where {T} = convert(T, 2G.n)
GroupsCore.gens(G::DihedralGroup) = [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]

# Base.rand not needed for our purposes here

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

# optional functions:
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

# The Robinson form is invariant under the following action of the Dihedral group on monomials:
# The action of each element of the groups is to map the variables `x, y` to:
#
# | id | rotation | reflection |
# |----|----------|------------|
# | 0  | x, y     | y, x       |
# | 1  | -y, x    | -x, y      |
# | 2  | -x, -y   | -y, -x     |
# | 3  | y, -x    | x, -y      |

using SumOfSquares
using DynamicPolynomials
@polyvar x y
struct DihedralAction <: Symmetry.OnMonomials end
import SymbolicWedderburn
SymbolicWedderburn._coeff_type(::DihedralAction) = Float64
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

# We can verify that `poly` is indeed invariant under the action of each element of the group as follows.

G = DihedralGroup(4)
for g in G
    @show SymbolicWedderburn.action(DihedralAction(), g, poly)
    @test SymbolicWedderburn.action(DihedralAction(), g, poly) == poly #src
end

# We can exploit this symmetry for reducing the problem using the `SymmetricIdeal` certificate as follows:

import CSDP
function solve(G)
    solver = CSDP.Optimizer
    model = Model(solver)
    @variable(model, t)
    @objective(model, Max, t)
    pattern = Symmetry.Pattern(G, DihedralAction())
    con_ref = @constraint(model, poly - t in SOSCone(), symmetry = pattern)
    optimize!(model)
    @test value(t) ≈ -3825/4096 rtol=1e-2 #src
    @show value(t)


    g = gram_matrix(con_ref).sub_gram_matrices #src
    @test length(g) == 5                       #src
    @test g[1].basis.polynomials == [x*y^2, x, x^3] #src
    @test g[2].basis.polynomials == [x^2*y, y, y^3] #src
    for i in 1:2                        #src
        I = 3:-1:1                      #src
        Q = g[i].Q[I, I]                #src
        @test size(Q) == (3, 3)         #src
        @test Q[1, 1] ≈  1    rtol=1e-3 #src
        @test Q[1, 2] ≈ -5/8  rtol=1e-3 #src
        @test Q[1, 3] ≈ -1    rtol=1e-3 #src
        @test Q[2, 2] ≈ 25/64 rtol=1e-3 #src
        @test Q[2, 3] ≈  5/8  rtol=1e-3 #src
        @test Q[3, 3] ≈  1    rtol=1e-3 #src
    end #src
    @test g[3].basis.polynomials == [x^2 + y^2, 1.0] #src
    @test size(g[3].Q) == (2, 2)             #src
    @test g[3].Q[2, 2] ≈ 7921/4096 rtol=1e-3 #src
    @test g[3].Q[1, 2] ≈  -89/128 rtol=1e-3  #src
    @test g[3].Q[1, 1] ≈ 1/4 rtol=1e-2 #src
    @test g[4].basis.polynomials == [x * y] #src
    @test size(g[4].Q) == (1, 1)       #src
    @test g[4].Q[1, 1] ≈ 0   atol=1e-3 #src
    @test g[5].basis.polynomials == [x^2 - y^2] #src
    @test size(g[5].Q) == (1, 1)       #src
    @test g[5].Q[1, 1] ≈ 0   atol=1e-3 #src
    for g in gram_matrix(con_ref).sub_gram_matrices
        println(g.basis.polynomials)
    end
end
solve(G)

# We notice that we indeed find `-3825/4096` and that symmetry was exploited.
# In case the conjugacy classes are known, we can implement
# `SymbolicWedderburn.conjugacy_classes_orbit` instead of `order` and `iterate`.
# To show that these do not need to be implemented, we create a new dihedral group type
# that do not implement these methods but that instead implement
# `SymbolicWedderburn.conjugacy_classes_orbit`:

struct DihedralGroup2 <: GroupsCore.Group
    n::Int
end
PermutationGroups.gens(G::DihedralGroup2) = [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]
_orbit(cc::Vector{<:GroupsCore.GroupElement}) = PermutationGroups.Orbit(cc, Dict(a => nothing for a in cc))
_orbit(el::GroupsCore.GroupElement) = _orbit([el])
function SymbolicWedderburn.conjugacy_classes_orbit(d::DihedralGroup2)
    orbits = [_orbit(DihedralElement(d.n, false, 0))]
    for i in 1:div(d.n - 1, 2)
        push!(orbits, _orbit([
            DihedralElement(d.n, false, i),
            DihedralElement(d.n, false, d.n - i),
        ]))
    end
    if iseven(d.n)
        push!(orbits, _orbit(DihedralElement(d.n, false, div(d.n, 2))))
        push!(orbits, _orbit([
            DihedralElement(d.n, true, i) for i in 0:2:(d.n - 2)
        ]))
        push!(orbits, _orbit([
            DihedralElement(d.n, true, i) for i in 1:2:(d.n - 1)
        ]))
    else
        push!(orbits, _orbit([
            DihedralElement(d.n, true, i) for i in 0:(d.n - 1)
        ]))
    end
end

# As we have not implemented the iterator over all the elements, we can iterate over all
# conjugacy classes instead to verify that the polynomial is invariant under the group action.

G = DihedralGroup2(4)
for cc in SymbolicWedderburn.conjugacy_classes_orbit(G)
    for g in cc
        @show SymbolicWedderburn.action(DihedralAction(), g, poly)
        @test SymbolicWedderburn.action(DihedralAction(), g, poly) == poly #src
    end
end

solve(G)
