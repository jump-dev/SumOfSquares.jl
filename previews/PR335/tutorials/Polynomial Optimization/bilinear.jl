# # Bilinear terms

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/bilinear.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/bilinear.ipynb)
# **Adapted from**: Section 3.1 of [F99] and Table 5.1 of [Las09]
#
# [F99] Floudas, Christodoulos A. et al.
# *Handbook of Test Problems in Local and Global Optimization.*
# Nonconvex Optimization and Its Applications (NOIA, volume 33).
#
# [Las09] Lasserre, J. B.
# *Moments, positive polynomials and their applications*
# World Scientific, **2009**.

# ## Introduction

# Consider the polynomial optimization problem [F99, Section 3.1]

using Test #src
using DynamicPolynomials
@polyvar x[1:8]
p = sum(x[1:3])
using SumOfSquares
K = @set 0.0025 * (x[4] + x[6]) <= 1 &&
    0.0025 * (-x[4] + x[5] + x[7]) <= 1 &&
    0.01 * (-x[5] + x[8]) <= 1 &&
    100x[1] - x[1] * x[6] + 8333.33252x[4] <= 250000/3 &&
    x[2] * x[4] - x[2] * x[7] - 1250x[4] + 1250x[5] <= 0 &&
    x[3] * x[5] - x[3] * x[8] - 2500x[5] + 1250000 <= 0 &&
    100 <= x[1] && x[1] <= 10000 &&
    1000 <= x[2] && x[2] <= 10000 &&
    1000 <= x[3] && x[3] <= 10000 &&
    10 <= x[4] && x[4] <= 1000 &&
    10 <= x[5] && x[5] <= 1000 &&
    10 <= x[6] && x[6] <= 1000 &&
    10 <= x[7] && x[7] <= 1000 &&
    10 <= x[8] && x[8] <= 1000

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

# The first level of the hierarchy gives a lower bound of `2100`

model2 = solve(2)
nothing # hide
@test objective_value(model2) ≈ 2100 rtol=1e-4 #src
@test termination_status(model2) == MOI.OPTIMAL #src
