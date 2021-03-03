# # Dihedral symmetry of the Robinson form

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Dihedral symmetry of the Robinson form.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Dihedral symmetry of the Robinson form.ipynb)
# **Adapted from**: Example 5.4 of [GP04]
#
# [GP04] Gatermann, Karin and Parrilo, Pablo A.
# *Symmetry groups, semidefinite programs, and sums of squares*.
# Journal of Pure and Applied Algebra 192.1-3 (2004): 95-128.


using Test #src

# Symmetry reduction is still a work in progress in SumOfSquares, so we include the following files that will be incorporated into SumOfSquares.jl once SymbolicWedderburn.jl is released:
using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "symmetry.jl"))
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "scaled_perm.jl"))

# We start by defining the Dihedral group of order 8.
# This group is isomorphic to the following permutation group:

d = perm"(1, 2, 3, 4)"
c = perm"(1, 3)"
G = PermGroup([c, d])

# We could rely on this isomorphism to define this group.
# However, in order to illustrate how to do symmetry reduction with a custom group,
# we show in this example what should be implemented to define a new group.

struct DihedralElement <: GroupElem
    n::Int
    reflection::Bool
    id::Int
end
function PermutationGroups.order(el::DihedralElement)
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
Base.one(el::Union{DihedralElement}) = DihedralElement(el.n, false, 0)
function Base.inv(el::DihedralElement)
    if el.reflection || iszero(el.id)
        return el
    else
        return DihedralElement(el.n, false, el.n - el.id)
    end
end
function PermutationGroups.mul!(::DihedralElement, a::DihedralElement, b::DihedralElement)
    a.n == b.n || error("Cannot multiply elements from different Dihedral groups")
    id = mod(a.reflection ? a.id - b.id : a.id + b.id, a.n)
    return DihedralElement(a.n, a.reflection != b.reflection, id)
end
function Base.:^(el::DihedralElement, k::Integer)
    if el.reflection
        return iseven(k) ? one(el) : el
    else
        return DihedralElement(el.n, false, mod(el.id * k, el.n))
    end
end

struct DihedralGroup <: Group
    n::Int
end
_orbit(cc::Vector{<:GroupElem}) = PermutationGroups.Orbit(cc, Dict(a => nothing for a in cc))
_orbit(el::GroupElem) = _orbit([el])
function SymbolicWedderburn.conjugacy_classes_orbit(d::DihedralGroup)
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

# The Robinson form is invariant under the following action of the Dihedral group on monomials:
# The action of each element of the groups is to map the variables `x, y` to:
#
# | id | rotation | reflection |
# |----|----------|------------|
# | 0  | x, y     | y, x       |
# | 1  | -y, x    | -x, y      |
# | 2  | -x, -y   | -y, -x     |
# | 3  | y, -x    | x, -y      |

using DynamicPolynomials
@polyvar x y
function action(mono::MP.AbstractMonomial, el::DihedralElement)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return MP.substitute(MP.Eval(), mono, [x, y] => [sign_x * var_x, sign_y * var_y])
end
function action(term::MP.AbstractTerm, el::DihedralElement)
    return MP.coefficient(term) * action(MP.monomial(term), el)
end
function action(poly::MP.AbstractPolynomial, el::DihedralElement)
    return MP.polynomial([action(term, el) for term in MP.terms(poly)])
end

poly = x^6 + y^6 - x^4 * y^2 - y^4 * x^2 - x^4 - y^4 - x^2 - y^2 + 3x^2 * y^2 + 1

# We can verify that `poly` is indeed invariant under each element of the group as follows.

G = DihedralGroup(4)
for cc in SymbolicWedderburn.conjugacy_classes_orbit(G)
    for g in cc
        @show action(poly, g)
        @show action(poly, g) == poly
    end
end

# We can exploit this symmetry for reducing the problem using the `SymmetricIdeal` certificate as follows:

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G, action))
optimize!(model)
@test value(t) ≈ -3825/4096 rtol=1e-3 #src
value(t)

# We indeed find `-3825/4096`, let's verify that symmetry was exploited:

g = gram_matrix(con_ref).sub_gram_matrices #src
@test length(g) == 4                       #src
@test g[1].basis.polynomials == [x^3, x^2*y, x*y^2, y^3, x, y] #src
for I in [[5, 1, 3], [6, 4, 2]]            #src
    Q = g[1].Q[I, I]                       #src
    @test Q[1, 1] ≈ 25/64 rtol=1e-3 #src
    @test Q[1, 2] ≈ -5/8  rtol=1e-3 #src
    @test Q[1, 3] ≈  5/8  rtol=1e-3 #src
    @test Q[2, 2] ≈  1    rtol=1e-3 #src
    @test Q[2, 3] ≈ -1    rtol=1e-3 #src
    @test Q[3, 3] ≈  1    rtol=1e-3 #src
end #src
@test g[2].basis.polynomials == [x^2 + y^2, 1.0] #src
@test g[2].Q[2, 2] ≈ 7921/4096 rtol=1e-3 #src
@test g[2].Q[1, 2] ≈  -89/128 rtol=1e-3  #src
@test g[2].Q[1, 1] ≈ 1/4 rtol=1e-3 #src
@test g[3].basis.polynomials == [x * y] #src
@test g[3].Q[1, 1] ≈ 0   atol=1e-3 #src
@test g[4].basis.polynomials == [x^2 - y^2] #src
@test g[4].Q[1, 1] ≈ 0   atol=1e-3 #src
for g in gram_matrix(con_ref).sub_gram_matrices
    println(g.basis.polynomials)
end
