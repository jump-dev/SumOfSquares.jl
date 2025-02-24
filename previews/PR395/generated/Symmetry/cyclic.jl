using GroupsCore
import PermutationGroups

struct CyclicElem <: GroupElement
    n::Int
    id::Int
end
Base.:(==)(a::CyclicElem, b::CyclicElem) = a.n == b.n && a.id == b.id
Base.inv(el::CyclicElem) = CyclicElem(el.n, (el.n - el.id) % el.n)

function Base.:*(a::CyclicElem, b::CyclicElem)
    return CyclicElem(a.n, (a.id + b.id) % a.n)
end
Base.:^(el::CyclicElem, k::Integer) = CyclicElem(el.n, (el.id * k) % el.n)

Base.conj(a::CyclicElem, b::CyclicElem) = inv(b) * a * b
Base.:^(a::CyclicElem, b::CyclicElem) = conj(a, b)

function PermutationGroups.order(el::CyclicElem)
    return div(el.n, gcd(el.n, el.id))
end

struct CyclicGroup <: Group
    n::Int
end
Base.eltype(::CyclicGroup) = CyclicElem
Base.one(c::Union{CyclicGroup, CyclicElem}) = CyclicElem(c.n, 0)
PermutationGroups.gens(c::CyclicGroup) = [CyclicElem(c.n, 1)]
PermutationGroups.order(::Type{T}, c::CyclicGroup) where {T} = convert(T, c.n)
function Base.iterate(c::CyclicGroup, prev::CyclicElem=CyclicElem(c.n, -1))
    id = prev.id + 1
    if id >= c.n
        return nothing
    else
        next = CyclicElem(c.n, id)
        return next, next
    end
end

import MultivariatePolynomials as MP
import MultivariateBases as MB

using SumOfSquares

struct Action{V<:MP.AbstractVariable} <: Symmetry.OnMonomials
    variables::Vector{V}
end
Symmetry.SymbolicWedderburn.coeff_type(::Action) = Float64
function Symmetry.SymbolicWedderburn.action(a::Action, el::CyclicElem, mono::MP.AbstractMonomial)
    return prod(MP.powers(mono), init=MP.constant_monomial(mono)) do (var, exp)
        index = findfirst(isequal(var), a.variables)
        new_index = mod1(index + el.id, el.n)
        return a.variables[new_index]^exp
    end
end

using DynamicPolynomials
@polyvar x[1:3]
action = Action(x)
g = CyclicElem(3, 1)
Symmetry.SymbolicWedderburn.action(action, g, x[1]^3 * x[2] * x[3]^4)

N = 3
G = CyclicGroup(N)
poly = sum(x[i] * x[mod1(i + 1, N)] for i in 1:N) + sum(x.^2)
Symmetry.SymbolicWedderburn.action(action, g, poly)

import Clarabel
solver = Clarabel.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, action)
basis = MB.explicit_basis(MB.algebra_element(poly - t))
using SymbolicWedderburn
summands = SymbolicWedderburn.symmetry_adapted_basis(
    Float64,
    pattern.group,
    pattern.action,
    basis,
    semisimple = true,
)

gram_basis = SumOfSquares.Certificate.Symmetry._gram_basis(pattern, basis, Float64)

con_ref = @constraint(model, poly - t in SOSCone(), symmetry = pattern)
optimize!(model)
solution_summary(model)

gram_matrix(con_ref)

basis = [(x[1] + x[2] - 2x[3])/√6, (x[1] - x[2])/√2]

image = [Symmetry.SymbolicWedderburn.action(action, g, p) for p in basis]

a = -1/2
b = √3/2
[a -b; b a] * basis

model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, action)
cone = SumOfSquares.NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
con_ref = @constraint(model, poly - t in cone, symmetry = pattern)
optimize!(model)
solution_summary(model)

gram_matrix(con_ref)

complex_basis = basis[1] + im * basis[2]
image = Symmetry.SymbolicWedderburn.action(action, g, complex_basis)

(a + b * im) * complex_basis

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
