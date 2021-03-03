using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "symmetry.jl"))
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "scaled_perm.jl"))

d = perm"(1, 2, 3, 4)"
c = perm"(1, 3)"
G = PermGroup([c, d])

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

G = DihedralGroup(4)
for cc in SymbolicWedderburn.conjugacy_classes_orbit(G)
    for g in cc
        @show action(poly, g)
        @show action(poly, g) == poly
    end
end

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G, action))
optimize!(model)
value(t)

for g in gram_matrix(con_ref).sub_gram_matrices
    println(g.basis.polynomials)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

