# # Even reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Even reduction.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Even reduction.ipynb)

using Test #src
using DynamicPolynomials
@polyvar x

# We would like to find the minimum value of the polynomial

poly = x^4 - 2x^2

using SumOfSquares

# Symmetry reduction is still a work in progress in SumOfSquares, so we include the following files that will be incorporated into SumOfSquares.jl once SymbolicWedderburn.jl is released:
using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "symmetry.jl"))
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "scaled_perm.jl"))

# We define the custom action as follows:

using PermutationGroups
function action(mono::MP.AbstractMonomial, p::Perm)
    if p == perm"(1)(2)" || iseven(MP.degree(mono))
        return 1 * mono
    else
        @assert p == perm"(1,2)"
        return -1 * mono
    end
end
G = PermGroup([perm"(1,2)"])

# We can exploit the symmetry as follows:

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

@test length(gram_matrix(con_ref).sub_gram_matrices) == 2 #src
@test gram_matrix(con_ref).sub_gram_matrices[1].basis.polynomials == [x^2, 1] #src
@test gram_matrix(con_ref).sub_gram_matrices[2].basis.polynomials == [x] #src
gram_matrix(con_ref)
