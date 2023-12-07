# # Nonconvex QP

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/qp.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/qp.ipynb)
# **Adapted from**: Section 2.2 of [F99], [Lasserre2009; Table 5.1](@cite)
#
# [F99] Floudas, Christodoulos A. et al.
# *Handbook of Test Problems in Local and Global Optimization.*
# Nonconvex Optimization and Its Applications (NOIA, volume 33).

# ## Introduction

# Consider the nonconvex Quadratic Program (QP) [F99, Section 2.2]
# that minimizes the *concave* function $c^\top x - x^\top Qx / 2$
# over the polyhedron obtained by intersecting the hypercube $[0, 1]^5$
# with the halfspace $10x_1 + 12x_2 + 11x_3 + 7x_4 + 4x_5 \le 40$.

using Test #src

using LinearAlgebra
c = [42, 44, 45, 47, 47.5]
Q = 100I

using DynamicPolynomials
@polyvar x[1:5]
p = c'x - x' * Q * x / 2
using SumOfSquares
K = @set x[1] >= 0 && x[1] <= 1 &&
         x[2] >= 0 && x[2] <= 1 &&
         x[3] >= 0 && x[3] <= 1 &&
         x[4] >= 0 && x[4] <= 1 &&
         x[5] >= 0 && x[5] <= 1 &&
         10x[1] + 12x[2] + 11x[3] + 7x[4] + 4x[5] <= 40

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.

import Clarabel
solver = Clarabel.Optimizer

# A Sum-of-Squares certificate that $p \ge \alpha$ over the domain `S`, ensures that $\alpha$ is a lower bound to the polynomial optimization problem.
# The following function searches for the largest lower bound and finds zero using the `d`th level of the hierarchy`.

function solve(d)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = K, maxdegree = d)
    optimize!(model)
    println(solution_summary(model))
    return model
end

# The first level of the hierarchy does not give any lower bound:

model2 = solve(2)
nothing # hide
@test termination_status(model2) == MOI.INFEASIBLE #src

# Indeed, as the constraints have degree 1 and their multipliers are SOS
# so they have an even degree, with `maxdegree` 2 we can only use degree 0
# multipliers hence constants. The terms of maximal degree in resulting
# sum will therefore only be in `-x' * Q * x/2` hence it is not SOS whatever
# is the value of the multipliers. Let's try with `maxdegree` 3 so that the
# multipliers can be quadratic.
# This second level is now feasible and gives a lower bound of `-22`.

model3 = solve(3)
nothing # hide
@test objective_value(model3) ≈ -22 rtol=1e-4 #src
@test termination_status(model3) == MOI.OPTIMAL #src
