# # Nonnegative over a variety

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/circle.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/circle.ipynb)

# The polynomial ``1 - y^2`` is nonnegative for all ``y`` in the unit circle.
# This can be verified using Sum-of-Squares.

using Test #src
using DynamicPolynomials
using SumOfSquares
@polyvar x y
S = @set x^2 + y^2 == 1

# We need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.
# The domain over which the nonnegativity of ``1 - y^2`` should be certified
# is specified through the `domain` keyword argument.

import CSDP
model = SOSModel(CSDP.Optimizer)
set_silent(model)
con_ref = @constraint(model, 1 - y^2 >= 0, domain = S)
optimize!(model)

# We can see that the model was feasible:

@test JuMP.termination_status(model) == MOI.OPTIMAL #src
@test JuMP.primal_status(model) == MOI.FEASIBLE_POINT #src
solution_summary(model)

# The certificate can be obtained as follows:

dec = sos_decomposition(con_ref, 1e-6) #src
@test length(dec) == 1 #src
@test first(dec) ≈ x rtol = 1e-6 #src
sos_decomposition(con_ref, 1e-6)

# It returns ``x^2`` which is a valid certificate as:
# $$ 1 - y^2 \equiv x^2 \pmod{x^2 + y^2 - 1} $$
