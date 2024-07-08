# # Cyclic symmetry for Sums of Hermitian Squares

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry/cyclic.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry/cyclic.ipynb)

using Test #src

# We start by defining the Cyclic group.

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

# Now we define that the cyclic group acts on monomial by permuting variables
# cyclically. So for instance, `CyclicElem(3, 1)` would transform
# `x_1^3*x_2*x_3^4` into `x_1^4*x_2^3*x_3`.

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
@test Symmetry.SymbolicWedderburn.action(action, g, x[1]^3 * x[2] * x[3]^4) == x[1]^4 * x[2]^3 * x[3] #src
Symmetry.SymbolicWedderburn.action(action, g, x[1]^3 * x[2] * x[3]^4)

# The following polynomial `poly` is invariant under the action of the group `G`.

N = 3
G = CyclicGroup(N)
poly = sum(x[i] * x[mod1(i + 1, N)] for i in 1:N) + sum(x.^2)
@test Symmetry.SymbolicWedderburn.action(action, g, poly) == poly #src
Symmetry.SymbolicWedderburn.action(action, g, poly)

# Let's now find the minimum of `p` by exploiting this symmetry.

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, action)
#import MultivariateBases as MB
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
@test termination_status(model) == MOI.OPTIMAL #src
@test objective_value(model) ≈ 0.0 atol=1e-4 #src
solution_summary(model)

# Let's look at the symmetry adapted basis used.

gram = gram_matrix(con_ref).blocks #src
@test length(gram) == 3                       #src
@test gram[1].Q ≈ [0 0; 0 2] #src
@test length(gram[1].basis.polynomials) == 2 #src
@test gram[1].basis.polynomials[1] == 1 #src
@test gram[1].basis.polynomials[2] ≈ -sum(x)/√3 #src
@test gram[2].Q ≈ [0.5;;] #src
@test length(gram[2].basis.polynomials) == 1 #src
@test gram[2].basis.polynomials[1] ≈ (x[1] + x[2] - 2x[3])/√6 #src
@test gram[3].Q == gram[2].Q #src
@test length(gram[3].basis.polynomials) == 1 #src
@test gram[3].basis.polynomials[1] ≈ (x[1] - x[2])/√2 #src
for gram in gram_matrix(con_ref).blocks
    println(gram.basis.polynomials)
    display(gram.Q)
end

# Let's look into more details at the last two elements of the basis.

basis = [(x[1] + x[2] - 2x[3])/√6, (x[1] - x[2])/√2]

# This actually constitutes the basis for an invariant subspace corresponding
# to a group character of degree 2 and multiplicity 1.
# This means that it decomposes the semidefinite matrix into 2 blocks of size
# 1-by-1 that are equal. Indeed, we see above that `gram.Q` is identically equal for both.
# As the group is generated by one element `g`, we can just verify it by verifying
# its invariance under `g`.
# The image of each element under the basis is:

image = [Symmetry.SymbolicWedderburn.action(action, g, p) for p in basis]

# We can see that they are both still in the same 2-dimensional subspace.

a = -1/2
b = √3/2
@test all(image .≈ [a -b; b a] * basis) #src
[a -b; b a] * basis

# In fact, these last two basis comes from the real decomposition of a complex one.

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, action)
cone = SumOfSquares.NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
con_ref = @constraint(model, poly - t in cone, symmetry = pattern)
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL #src
@test objective_value(model) ≈ 0.0 atol=1e-4 #src
solution_summary(model)

gram = gram_matrix(con_ref).blocks #src
@test length(gram) == 3                       #src
@test gram[1].Q ≈ [0 0; 0 2] #src
@test length(gram[1].basis.polynomials) == 2 #src
@test gram[1].basis.polynomials[1] == 1 #src
@test gram[1].basis.polynomials[2] ≈ -sum(x)/√3 #src
@test gram[2].Q ≈ [0.5;;] rtol = 1e-6 #src
@test length(gram[2].basis.polynomials) == 1 #src
@test gram[2].basis.polynomials[1] ≈ (basis[1] - basis[2] * im) / √2  #src
@test gram[3].Q ≈ [0.5;;] rtol = 1e-6 #src
@test length(gram[3].basis.polynomials) == 1 #src
@test gram[3].basis.polynomials[1] ≈ (basis[1] + basis[2] * im) / √2 #src
for gram in gram_matrix(con_ref).blocks
    println(gram.basis.polynomials)
    display(gram.Q)
end

# We can see that the real invariant subspace was in fact coming from two complex conjugate complex invariant subspaces:

complex_basis = basis[1] + im * basis[2]
image = Symmetry.SymbolicWedderburn.action(action, g, complex_basis)

# And there is a direct correspondance between the representation of the real and complex versions:

@test all(image .≈ (a + b * im) * complex_basis) #src
(a + b * im) * complex_basis
