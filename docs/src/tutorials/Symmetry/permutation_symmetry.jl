# # Symmetry reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry/symmetry_reduction.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry/symmetry_reduction.ipynb)
# **Adapted from**: [SymbolicWedderburn example](https://github.com/kalmarek/SymbolicWedderburn.jl/blob/tw/ex_sos/examples/ex_C4.jl)

import MutableArithmetics as MA
using MultivariatePolynomials
using MultivariateBases

using Test #src
using DynamicPolynomials
@polyvar x[1:4]

# We would like to find the minimum value of the polynomial

poly = sum(x) + sum(x.^2)

# As we can decouple the problem for each `x[i]` for which `x[i] + x[i]^2` has
# minimum value 0.25, we would expect to get `-1` as answer.
# Can this decoupling be exploited by SumOfSquares as well ?
# For this, we need to use a certificate that can exploit the permutation symmetry of the polynomial.

using SumOfSquares

# We define the symmetry group as a permutation group in the variables.
# In order to do that, we define the action of a permutation on a monomial
# as the monomial obtained after permuting the variables.

using PermutationGroups
G = PermGroup([perm"(1,2,3,4)"])

# We can use this certificate as follows:

import Clarabel
solver = Clarabel.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, Symmetry.VariablePermutation())
con_ref = @constraint(model, poly - t in SOSCone(), symmetry = pattern)
optimize!(model)
@test value(t) ≈ -1 rtol=1e-6 #src
value(t)

# We indeed find `-1`, let's verify that symmetry was exploited:

gram = gram_matrix(con_ref).blocks    #src
@test length(gram) == 3                          #src
# `gram_basis` now returns one `SimpleBasis` per character (was `SemisimpleBasis`). #src
# Q sizes still match the multiplicity of each irreducible character. #src
@test size(gram[1].Q) == (2, 2)                  #src
@test size(gram[2].Q) == (1, 1)                  #src
@test size(gram[3].Q) == (1, 1)                  #src
gram_matrix(con_ref)
