# # Exterior of ellipsoid

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/ellipsoid.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/ellipsoid.ipynb)
# **Adapted from**: Section 3.5 of [F99], Example 6.23 of [L09] and Table 5.1 of [Las09]
#
# [F99] Floudas, Christodoulos A. et al.
# *Handbook of Test Problems in Local and Global Optimization.*
# Nonconvex Optimization and Its Applications (NOIA, volume 33).
#
# [Las09] Lasserre, J. B.
# *Moments, positive polynomials and their applications*
# World Scientific, **2009**.

# ## Introduction

# Consider the polynomial optimization problem [F99, Section 3.5]

A = [
     0  0  1
     0 -1  0
    -2  1 -1
]
bz = [3, 0, -4] - [0, -1, -6]
y = [1.5, -0.5, -5]

using Test #src
using DynamicPolynomials
@polyvar x[1:3]
p = -2x[1] + x[2] - x[3]
using SumOfSquares
e = A * x - y
f = e'e - bz'bz / 4
K = @set sum(x) <= 4 && 3x[2] + x[3] <= 6 && f >= 0 && 0 <= x[1] && x[1] <= 2 && 0 <= x[2] && 0 <= x[3] && x[3] <= 3

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

# The first level of the hierarchy gives a lower bound of `-7``

model2 = solve(2)
@test objective_value(model2) ≈ -6 rtol=1e-4 #src
@test termination_status(model2) == MOI.OPTIMAL #src

# The second level improves the lower bound

model3 = solve(4)
@test objective_value(model5) ≈ -74/13 rtol=1e-4 #src
@test termination_status(model5) == MOI.OPTIMAL #src

# The third level improves it even further

model6 = solve(6)
@test objective_value(model6) ≈ -4.06848 rtol=1e-4 #src
@test termination_status(model6) == MOI.OPTIMAL #src

# The fourth level finds the optimal objective value as lower bound.

model8 = solve(8)
@test objective_value(model8) ≈ -4 rtol=1e-4 #src
@test termination_status(model8) == MOI.OPTIMAL #src
