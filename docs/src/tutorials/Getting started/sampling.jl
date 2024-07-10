# # Sampling basis

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/sampling.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/sampling.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
using DynamicPolynomials
using SumOfSquares
import Hypatia

# In this tutorial, we show how to use a different polynomial basis
# for enforcing the equality between the polynomial and its Sum-of-Squares decomposition.

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 9

model = Model(Hypatia.Optimizer)
set_silent(model)
@constraint(model, p in SOSCone(), zero_basis = BoxSampling([-1], [1]))
optimize!(model)
solution_summary(model)
backend(model)
