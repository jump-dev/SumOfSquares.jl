# # Symmetry reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry reduction.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry reduction.ipynb)
# **Adapted from**: https://github.com/kalmarek/SymbolicWedderburn.jl/blob/tw/ex_sos/examples/ex_C4.jl

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using Test #src
using DynamicPolynomials
@polyvar x[1:4]

# We would like to find the minimum value of the polynomial

poly = sum(x) + sum(x.^2)

# As we can decouple the problem for each `x[i]` for which `x[i] + x[i]^2` has
# minimum value 0.25, we would expect to get `-1` as answer.
# Can this decoupling be exploited by SumOfSquares as well ?
# For this, we need to use a certificate that can exploit the permutation symmetry of the polynomial.
# Symmetry reduction is still a work in progress in SumOfSquares, so we include the following files that will be incorporated into SumOfSquares.jl once SymbolicWedderburn.jl is released:
using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "symmetry.jl"))

# We define the symmetry group as a permutation group in the variables.
# In order to do that, we define the action of a permutation on a monomial
# as the monomial obtained after permuting the variables.

using PermutationGroups
function action(mono::AbstractMonomial, p::Perm)
    v = variables(mono)
    MP.substitute(MP.Eval(), mono, v => [v[i^p] for i in eachindex(v)])
end
G = PermGroup([perm"(1,2,3,4)"])

# We can use this certificate as follows:

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G, action))
optimize!(model)
@test value(t) ≈ -1 #src
value(t)

# We indeed find `-1`, let's verify that symmetry was exploited:

@test length(gram_matrix(con_ref).sub_gram_matrices) == 3 #src
gram_matrix(con_ref)
