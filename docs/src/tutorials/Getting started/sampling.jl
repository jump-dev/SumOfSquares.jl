# # Sampling basis

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/sampling.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/sampling.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
using DynamicPolynomials
using SumOfSquares
import MultivariateBases as MB

using Dualization
import Hypatia
import SCS

# In this tutorial, we show how to use a sampling basis for enforcing
# the equality between the polynomial and its Sum-of-Squares decomposition.
# As introduced in [Zero basis](@ref), we can set the zero basis to a sampling
# basis using the `zero_basis` keyword.

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 3

model = Model(dual_optimizer(SCS.Optimizer))
set_silent(model)
@variable(model, γ)
@objective(model, Max, γ)
@constraint(model, p - γ in SOSCone(), zero_basis = BoxSampling([-1], [1]))
optimize!(model)
solution_summary(model)
@test primal_status(model) == MOI.FEASIBLE_POINT #src
@test value(γ) ≈ -6 rtol=1e-4 #src

# We can see that the SOS constraint is converted into a PSD constraint.

print_active_bridges(model)

# Let's try with Hypatia now:

set_optimizer(model, dual_optimizer(Hypatia.Optimizer))
optimize!(model)
solution_summary(model)
@test primal_status(model) == MOI.FEASIBLE_POINT #src
@test value(γ) ≈ -6 rtol=1e-4 #src

# We can see that the SOS constraint is passed as a
# `LowRankOpt.SetDotProducts{LowRankOpt.WITHOUT_SET}`. This maps into
# Hypatia's native `WSOSInterpNonnegativeCone`.

print_active_bridges(model)
