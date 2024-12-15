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
polys = gram[1].basis.bases[].elements #src
@test length(polys) == 2     #src
@test polys[1] ≈ 1        #src
@test polys[2] ≈ -0.5 * sum(x) #src
@test size(gram[1].Q) == (2, 2)                  #src
@test gram[1].Q[1, 1] ≈ 1.0 atol=1e-6            #src
@test gram[1].Q[1, 2] ≈ -1.0 atol=1e-6            #src
@test gram[1].Q[2, 2] ≈ 1.0 atol=1e-6            #src
@test length(gram[2].basis.bases) == 2 #src
polys = gram[2].basis.bases[1].elements #src
@test length(polys) == 1     #src
@test polys[] ≈ (x[2] - x[4]) / √2 #src
@test size(g[2].Q) == (1, 1)                  #src
@test g[2].Q[1, 1] ≈ 1.0 atol=1e-6            #src
polys = gram[2].basis.bases[2].elements #src
@test length(polys) == 1     #src
@test polys[1] ≈ (x[1] - x[3]) / √2 #src
polys = gram[3].basis.bases[].elements #src
@test length(polys) == 1     #src
@test polys[] ≈ (x[1] - x[2] + x[3] - x[4]) / 2 #src
@test size(gram[3].Q) == (1, 1)                  #src
@test gram[3].Q[1, 1] ≈ 1.0 atol=1e-6            #src
gram_matrix(con_ref)
