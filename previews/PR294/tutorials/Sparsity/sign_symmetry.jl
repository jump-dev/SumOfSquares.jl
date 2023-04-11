# # Sign symmetry

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Sparsity/sign_symmetry.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Sparsity/sign_symmetry.ipynb)
# **Adapted from**: Example 4 of [L09]
#
# [L09] Lofberg, Johan.
# *Pre-and post-processing sum-of-squares programs in practice*.
# IEEE transactions on automatic control 54, no. 5 (2009): 1007-1011.

using Test #src
using DynamicPolynomials
@polyvar x[1:3]

# We would like to determine whether the following polynomial is a sum-of-squares.

poly = 1 + x[1]^4 + x[1] * x[2] + x[2]^4 + x[3]^2

# In order to do this, we can solve the following Sum-of-Squares program.

import CSDP
solver = CSDP.Optimizer
using SumOfSquares
function sos_check(sparsity)
    model = Model(solver)
    con_ref = @constraint(model, poly in SOSCone(), sparsity = sparsity)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL #src
    println(solution_summary(model))
    return gram_matrix(con_ref)
end

g = sos_check(Sparsity.NoPattern())
@test g.basis.monomials == [x[1]^2, x[1] * x[2], x[2]^2, x[1], x[2], x[3], 1] #src
g.basis.monomials

# As detailed in the Example 4 of [L09], we can exploit the *sign symmetry* of
# the polynomial to decompose the large positive semidefinite matrix into smaller ones.

g = sos_check(Sparsity.SignSymmetry())
monos = [sub.basis.monomials for sub in g.sub_gram_matrices]
@test length(monos) == 3 #src
@test [x[1], x[2]] in monos #src
@test [x[3]] in monos #src
@test [x[1]^2, x[1] * x[2], x[2]^2, 1] in monos #src
