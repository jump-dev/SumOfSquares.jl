using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "symmetry.jl"))
#include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "scaled_perm.jl"))

d = perm"(1, 2, 3, 4)"
c = perm"(1, 3)"
G = PermGroup([c, d])

struct DihedralElement <: GroupElement
    n::Int
    reflection::Bool
    id::Int
end
function Base.:(==)(g::DihedralElement, h::DihedralElement)
    return g.n == h.n && g.reflection == h.reflection && g.id == h.id
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
Base.one(el::DihedralElement) = DihedralElement(el.n, false, 0)
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
function Base.:^(el::DihedralElement, k::Integer)
    if el.reflection
        return iseven(k) ? one(el) : el
    else
        return DihedralElement(el.n, false, mod(el.id * k, el.n))
    end
end

Base.conj(a::DihedralElement, b::DihedralElement) = inv(b) * a * b
Base.:^(a::DihedralElement, b::DihedralElement) = conj(a, b)

struct DihedralGroup <: Group
    n::Int
end
Base.one(G::DihedralGroup) = DihedralElement(G.n, false, 0)
PermutationGroups.gens(G::DihedralGroup) = [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]
PermutationGroups.order(::Type{T}, G::DihedralGroup) where {T} = convert(T, 2G.n)
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

using DynamicPolynomials
@polyvar x y
struct DihedralAction <: SymbolicWedderburn.ByLinearTransformation end
SymbolicWedderburn.coeff_type(::DihedralAction) = Float64
function SymbolicWedderburn.action(::DihedralAction, el::DihedralElement, mono::MP.AbstractMonomial)
    if iseven(el.reflection + el.id)
        var_x, var_y = x, y
    else
        var_x, var_y = y, x
    end
    sign_x = 1 <= el.id <= 2 ? -1 : 1
    sign_y = 2 <= el.id ? -1 : 1
    return MP.substitute(MP.Eval(), mono, [x, y] => [sign_x * var_x, sign_y * var_y])
end
function SymbolicWedderburn.action(a::DihedralAction, el::DihedralElement, term::MP.AbstractTerm)
    return MP.coefficient(term) * SymbolicWedderburn.action(a, el, MP.monomial(term))
end
function SymbolicWedderburn.action(a::DihedralAction, el::DihedralElement, poly::MP.AbstractPolynomial)
    return MP.polynomial([SymbolicWedderburn.action(a, el, term) for term in MP.terms(poly)])
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
    certificate = SymmetricIdeal(Certificate.MaxDegree(SOSCone(), MonomialBasis, maxdegree(poly)), G, DihedralAction())
    con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = certificate)
    optimize!(model)
    @show value(t)


    for g in gram_matrix(con_ref).sub_gram_matrices
        println(g.basis.polynomials)
    end
end
solve(G)

struct DihedralGroup2 <: Group
    n::Int
end
PermutationGroups.gens(G::DihedralGroup2) = [DihedralElement(G.n, false, 1), DihedralElement(G.n, true, 0)]
_orbit(cc::Vector{<:GroupElement}) = PermutationGroups.Orbit(cc, Dict(a => nothing for a in cc))
_orbit(el::GroupElement) = _orbit([el])
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

G = DihedralGroup2(4)
for cc in SymbolicWedderburn.conjugacy_classes_orbit(G)
    for g in cc
        @show SymbolicWedderburn.action(DihedralAction(), g, poly)
    end
end

solve(G)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

