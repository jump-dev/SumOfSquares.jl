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

# To use Sum-of-Squares Programming, we first need to pick an SDP solver,
# see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.

import Clarabel
using Dualization
model = SOSModel(dual_optimizer(Clarabel.Optimizer))

# Create a JuMP decision variable for the lower bound

@variable(model, γ)

# `f(x)` is the Goldstein-Price function

f1 = x[1] + x[2] + 1
f2 = 19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2
f3 = 2*x[1] - 3*x[2]
f4 = 18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2
f = (1 + f1^2*f2) * (30 + f3^2*f4)

# Constraints `f(x) - γ` to be a sum of squares

con_ref = @constraint(model, f >= γ)
@objective(model, Max, γ)
optimize!(model)

# The lower bound found is 3

@test objective_value(model) ≈ 3 rtol=1e-3 #src
solution_summary(model)

# The moment matrix is as follows, we can already see the global minimizer
# `[0, -1]` from the entries `(2, 1)` and `(3, 1)`.
# This heuristic way to obtain solutions to the polynomial optimization problem
# is suggested in [Laurent2008, (6.15)](@cite).

ν = moment_matrix(con_ref)

# Many entries of the matrix actually have the same moment.
# We can obtain the following list of these moments without duplicates
# (ignoring when difference of entries representing the same moments is below `1e-5`)

μ = measure(ν, atol = 1e-5)

# The truncated moment matrix can then be obtained as follows

ν_truncated = moment_matrix(μ, monomials(x, 0:3))

# Let's check if the flatness property is satisfied.
# The rank of `ν_truncated` seems to be 1:

using LinearAlgebra
LinearAlgebra.svdvals(Matrix(ν_truncated.Q))
LinearAlgebra.rank(Matrix(ν_truncated.Q), rtol = 1e-3)
@test LinearAlgebra.rank(Matrix(ν_truncated.Q), rtol = 1e-3) == 1 #src
svdvals(Matrix(ν_truncated.Q))

# The rank of `ν` is clearly higher than 1, closer to 3:

@test 3 <= LinearAlgebra.rank(Matrix(ν.Q), rtol = 1e-3) <= 4 #src
svdvals(Matrix(ν.Q))

# Even if the flatness property is not satisfied, we can
# still try extracting the minimizer with a low rank decomposition of rank 3.
# We find the optimal solution again doing so:

atoms = atomic_measure(ν, FixedRank(3)) #src
@test length(atoms.atoms) == 1 #src
@test atoms.atoms[1].center ≈ [0, -1] rtol=1e-3 #src
atomic_measure(ν, FixedRank(3))
