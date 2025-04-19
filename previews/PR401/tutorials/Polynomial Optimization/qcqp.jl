# # Nonconvex quadratically constrained quadratic programs

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/nonconvex_qcqp.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/nonconvex_qcqp.ipynb)
# **Adapted from**: [Hesse1973](@cite), [Floudas1999; Section 3.4](@cite), [Laurent2008; Example 6.22](@cite) and [Lasserre2009; Table 5.1](@cite)

# We consider the nonconvex Quadratically Constrained Quadratic Programs (QCQP)
# introduced in [H73].
# Consider now the polynomial optimization problem [Laurent2008; Example 6.22](@cite) of
# maximizing the convex quadratic function
# (hence nonconvex since convex programs should either maximize concave functions
# or minimize convex functions)
# $25(x_1 - 2)^2 + (x_2 - 2)^2 + (x_3 - 1)^2 + (x_4 - 4)^2 + (x_5 - 1)^2 + (x_6 - 4)^2$
# over the basic semialgebraic set defined by the nonconvex quadratic inequalities
# $(x_3 - 3)^2 + x_4 \ge 4$,
# $(x_5 - 3)^2 + x_6 \ge 4$,
# and linear inequalities
# $x_1 - 3x_2 \le 2$,
# $-x_1 + x_2 \le 2$,
# $2 \le x_1 + x_2 \le 6$,
# $0 \le x_1, x_2$,
# $1 \le x_3 \le 5$,
# $0 \le x_4 \le 6$,
# $1 \le x_5 \le 5$,
# $0 \le x_6 \le 10$,
# $x_2 \le 4x_1^4 - 32x_1^3 + 88x_1^2 - 96x_1 + 36$ and the box constraints
# $0 \le x_1 \le 3$ and $0 \le x_2 \le 4$,

using Test #src
using DynamicPolynomials
@polyvar x[1:6]
centers = [2, 2, 1, 4, 1, 4]
weights = [25, 1, 1, 1, 1, 1]
p = -weights' * (x .- centers).^2
using SumOfSquares
K = @set x[1] >= 0 && x[2] >= 0 &&
    x[3] >= 1 && x[3] <= 5 &&
    x[4] >= 0 && x[4] <= 6 &&
    x[5] >= 1 && x[5] <= 5 &&
    x[6] >= 0 && x[6] <= 10 &&
    (x[3] - 3)^2 + x[4] >= 4 &&
    (x[5] - 3)^2 + x[6] >= 4 &&
    x[1] - 3x[2] <= 2 &&
    -x[1] + x[2] <= 2 &&
    x[1] + x[2] <= 6 &&
    x[1] + x[2] >= 2

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

# The first level of the hierarchy cannot find any lower bound.

model2 = solve(2)
nothing # hide
@test termination_status(model2) == MOI.INFEASIBLE #src

# The second level of the hierarchy finds the lower bound of `-310`.

model3 = solve(4)
nothing # hide
@test termination_status(model3) == MOI.OPTIMAL #src
@test objective_value(model3) ≈ -310 rtol=1e-4 #src
