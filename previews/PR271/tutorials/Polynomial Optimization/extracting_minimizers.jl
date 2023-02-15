# # Extracting minimizers

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/extracting_minimizers.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/extracting_minimizers.ipynb)
# **Adapted from**: Example 6.23 of [L09]
#
# [L09] Laurent, Monique.
# *Sums of squares, moment matrices and optimization over polynomials.*
# Emerging applications of algebraic geometry (2009): 157-270.

# ## Introduction

# Consider the polynomial optimization problem [L09, Example 6.23] of
# minimizing the linear function $-x_1 - x_2$
# over the basic semialgebraic set defined by the inequalities
# $x_2 \le 2x_1^4 - 8x_1^3 + 8x_1^2 + 2$,
# $x_2 \le 4x_1^4 - 32x_1^3 + 88x_1^2 - 96x_1 + 36$ and the box constraints
# $0 \le x_1 \le 3$ and $0 \le x_2 \le 4$,
# World Scientific, **2009**.

using Test #src
using DynamicPolynomials
@polyvar x[1:2]
p = -sum(x)
using SumOfSquares
K = @set x[1] >= 0 && x[1] <= 3 && x[2] >= 0 && x[2] <= 4 && x[2] <= 2x[1]^4 - 8x[1]^3 + 8x[1]^2 + 2 && x[2] <= 4x[1]^4 - 32x[1]^3 + 88x[1]^2 - 96x[1] + 36

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v0.21.6/installation/#Supported-solvers) for a list of the available choices.

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

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

model4 = solve(4)
@test objective_value(model4) ≈ -7 rtol=1e-4 #src
@test termination_status(model4) == MOI.OPTIMAL #src

# The second level improves the lower bound but still does not provide the solution

model5 = solve(5)
@test objective_value(model5) ≈ -20/3 rtol=1e-4 #src
@test termination_status(model5) == MOI.OPTIMAL #src
ν5 = moment_matrix(model5[:c])
@test extractatoms(ν5, 1e-2) === nothing #src
@test extractatoms(ν5, 1e-3) === nothing #src
@test extractatoms(ν5, 1e-4) === nothing #src
@test extractatoms(ν5, 1e-5) === nothing #src
extractatoms(ν5, 1e-2)
SumOfSquares.MultivariateMoments.computesupport!(ν5, 1e-1)
ν5.support

# The third level finds the optimal objective value as lower bound...

model7 = solve(7)
@test objective_value(model7) ≈ -5.5080 rtol=1e-4 #src
@test termination_status(model7) == MOI.OPTIMAL #src

# ...and proves it by exhibiting the minimizer.

ν7 = moment_matrix(model7[:c])
@test extractatoms(ν3, 1e-3) === nothing #src
η = extractatoms(ν7, 1e-3) # Returns nothing as the dual is not atomic
@test length(η.atoms) == 1 #src
@test η.atoms[1].center ≈ [2.3295, 3.1785] rtol=1e-4 #src

# We can indeed verify that the objective value at `x_opt` is equal to the lower bound.

x_opt = η.atoms[1].center
@test x_opt ≈ [2.3295, 3.1785] rtol=1e-4 #src
p(x_opt)
