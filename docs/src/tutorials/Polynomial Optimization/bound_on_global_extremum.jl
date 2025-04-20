# # Bound on Global Extremum

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/bound_on_blobal_extremum.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/bound_on_global_extremum.ipynb)
# **Adapted from**: SOSTOOLS' SOSDEMO3 (See Section 4.3 of [SOSTOOLS User's Manual](https://github.com/oxfordcontrol/SOSTOOLS/blob/SOSTOOLS400/docs/SOSTOOLS_400.pdf))

using Test #src
using DynamicPolynomials
@polyvar x1 x2

# The Goldstein-Price function $f(x)$ is defined as follows:

f1 = x1 + x2 + 1
f2 = 19 - 14x1 + 3x1^2 - 14x2 + 6x1*x2 + 3x2^2
f3 = 2x1 - 3x2
f4 = 18 - 32x1 + 12x1^2 + 48x2 - 36x1*x2 + 27x2^2
f = (1 + f1^2 * f2) * (30 + f3^2 * f4)

# We need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.
# We use `SOSModel` instead of `Model` to be able to use the `>=` syntax for Sum-of-Squares constraints.

using SumOfSquares
using CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

# We create the decision variable $\gamma$ that will be the lower bound to the Goldstein-Price function.
# We maximize it to have the highest possible lower bound.

@variable(model, γ)
@objective(model, Max, γ)

# We constrain $\gamma$ to be a lower bound with the following constraint
# that ensures that $f(x_1, x_2) \ge \gamma$ for all $x_1, x_2$.

@constraint(model, f >= γ)

JuMP.optimize!(model)

# We verify that the solver has found a feasible solution:

JuMP.primal_status(model)
@test JuMP.primal_status(model) == MOI.FEASIBLE_POINT || JuMP.primal_status(model) == MOI.NEARLY_FEASIBLE_POINT #src

# We can now obtain the lower bound either with `value(γ)` or `objective_value(model)`:

objective_value(model)
@test objective_value(model) ≈ 3 rtol=1e-2 #src
