# # Goldstein-price function

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/goldstein_price.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/goldstein_price.ipynb)
# **Contributed by**: Benoît Legat

# In this example, we consider the minimization of the [Goldstein-price function](https://en.wikipedia.org/wiki/Test_functions_for_optimization).

using Test #src
using SumOfSquares
using DynamicPolynomials

# Create *symbolic* variables (not JuMP *decision* variables)

@polyvar x[1:2]

# To use Sum of Squares Programming, we first need to pick an SDP solver,
# see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.

import SCS
scs = SCS.Optimizer
import Dualization
dual_scs = Dualization.dual_optimizer(scs)
model = SOSModel(dual_scs)

# Create a JuMP decision variable for the lower bound

@variable(model, γ)

# f(x) is the Goldstein-Price function

f1 = x[1] + x[2] + 1
f2 = 19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2
f3 = 2*x[1] - 3*x[2]
f4 = 18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2
f = (1 + f1^2*f2) * (30 + f3^2*f4)

# Constraints f(x) - γ to be sum of squares

con_ref = @constraint(model, f >= γ)
@objective(model, Max, γ)
optimize!(model)

# The lower bound found is 3

@test objective_value(model) ≈ 3 rtol=1e-4 #src
solution_summary(model)

# The moment matrix is as follows, we can already see the global minimizer
# `[0, -1]` from the entries `(2, 1)` and `(3, 1)`.
# This heuristic way to obtain solutions to the polynomial optimization problem
# is suggested in [L09, (6.15)].
#
# [L09] Laurent, Monique.
# *Sums of squares, moment matrices and optimization over polynomials.*
# Emerging applications of algebraic geometry (2009): 157-270.


M = moment_matrix(con_ref)
