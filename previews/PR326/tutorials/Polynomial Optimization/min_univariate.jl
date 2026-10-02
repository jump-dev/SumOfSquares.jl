# # Maximizing as minimum

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/min_univariate.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/min_univariate.ipynb)
# **Adapted from**: [Floudas1999; Section 4.10](@cite), [Laurent2008; Example 6.23](@cite) and [Lasserre2009; Table 5.1](@cite)

# ## Introduction

# Consider the polynomial optimization problem from [Floudas1999; Section 4.10](@cite)
# of minimizing the linear function $-x_1 - x_2$
# over the basic semialgebraic set defined by the inequalities
# $x_2 \le 2x_1^4 - 8x_1^3 + 8x_1^2 + 2$,
# $x_2 \le 4x_1^4 - 32x_1^3 + 88x_1^2 - 96x_1 + 36$ and the box constraints
# $0 \le x_1 \le 3$ and $0 \le x_2 \le 4$,

using Test #src
using DynamicPolynomials
@polyvar x[1:2]
p = -sum(x)
using SumOfSquares
f1 = 2x[1]^4 - 8x[1]^3 + 8x[1]^2 + 2
f2 = 4x[1]^4 - 32x[1]^3 + 88x[1]^2 - 96x[1] + 36
K = @set x[1] >= 0 && x[1] <= 3 && x[2] >= 0 && x[2] <= 4 && x[2] <= f1 && x[2] <= f2

# As we can observe below, the bounds on `x[2]` could be dropped and
# optimization problem is equivalent to the maximization of `min(f1, f2)`
# between `0` and `3`.

xs = range(0, stop = 3, length = 100)
using Plots
plot(xs, f1.(xs), label = "f1")
plot!(xs, f2.(xs), label = "f2")
plot!(xs, 4 * ones(length(xs)), label = nothing)

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) for a list of the available choices.

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

model4 = solve(4)
nothing # hide
@test objective_value(model4) ≈ -7 rtol=1e-4 #src
@test termination_status(model4) == MOI.OPTIMAL #src

# The moment matrix is not flat so we cannot extract a minimizer with
# `atomic_measure`. We can still look for a feasible solution close to the
# vector of first-order moments ``(\mathbb{E}[x_1], \mathbb{E}[x_2])``,
# as suggested in [Laurent2008](@cite).
# This vector is not necessarily feasible so [`round_solution`](@ref)
# projects it onto `K` with [`heuristic_projection`](@ref).
# The constraints `0 ≤ x[1] ≤ 3` only depend on `x[1]` so `x[1]` is first
# projected onto `[0, 3]`. Once `x[1]` is fixed, the remaining constraints only
# depend on `x[2]` so `x[2]` is then projected onto `[0, min(4, f1(x[1]), f2(x[1]))]`.

ν4 = moment_matrix(model4[:c])
x4 = round_solution(ν4, K, p)
@test x4 ≈ [3, 0] atol=1e-6 #src

# The objective value at this feasible point is an upper bound to the optimal
# objective value so we now know that it is between `-7` and `-3`.

p(x4)

# Instead of only considering the first-order moments, we can sample
# from the Gaussian distribution with the same moments of order up to 2,
# project the samples and keep the best feasible point found.
# This generalizes the random hyperplane rounding of [Goemans1995](@cite),
# see also [Barak2016](@cite). We fix the seed of the random number generator
# so that the results are reproducible.

import Random
gaussian = GaussianRounding(rng = Random.MersenneTwister(0))
p(round_solution(ν4, K, p, rounding = gaussian))

# The second level improves the lower bound

model5 = solve(5)
nothing # hide
@test objective_value(model5) ≈ -20/3 rtol=1e-4 #src
@test termination_status(model5) == MOI.OPTIMAL #src

# The upper bound obtained from the first-order moments is improved as well:

ν5 = moment_matrix(model5[:c])
x5 = round_solution(ν5, K, p)
@test p(x5) ≈ -3.9012 rtol=1e-3 #src
p(x5)

# With the Gaussian rounding, we obtain a better upper bound:

x5_gaussian = round_solution(ν5, K, p, rounding = gaussian)
@test objective_value(model5) <= p(x5_gaussian) <= p(x5) #src
p(x5_gaussian)

# The third level finds the optimal objective value as lower bound...

model7 = solve(7)
nothing # hide
@test objective_value(model7) ≈ -5.5080 rtol=1e-4 #src
@test termination_status(model7) == MOI.OPTIMAL #src

# ...and proves it by exhibiting the minimizer.

ν7 = moment_matrix(model7[:c])
η = atomic_measure(ν7, 1e-3)
@test length(η.atoms) == 1 #src
@test η.atoms[1].center ≈ [2.3295, 3.1785] rtol=1e-4 #src

# The first-order moments now coincide with the minimizer so
# [`round_solution`](@ref) finds it too:

@test round_solution(ν7, K, p) ≈ [2.3295, 3.1785] rtol=1e-4 #src
round_solution(ν7, K, p)

# We can indeed verify that the objective value at `x_opt` is equal to the lower bound.

x_opt = η.atoms[1].center
@test x_opt ≈ [2.3295, 3.1785] rtol=1e-4 #src
p(x_opt)

# We can see visualize the solution as follows:

scatter!([x_opt[1]], [x_opt[2]], markershape = :star, label = nothing)
