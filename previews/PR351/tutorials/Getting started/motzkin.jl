# # Motzkin

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/motzkin.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/motzkin.ipynb)
# **Adapted from**: [Blekherman2012; (3.6) and (3.19)](@cite)

# The first explicit example of nonnegative polynomial that is not a sum of squares was found by Motzkin in 1967. By the [Arithmetic-geometric mean](https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean),
# $$ \frac{x^4y^2 + x^2y^4 + 1}{3} \ge \sqrt[3]{x^4y^2 \cdot x^2y^4 \cdot 1} = x^2y^2 $$
# hence
# $$ x^4y^2 + x^2y^4 + 1 - 3x^2y^2 \ge 0. $$
# The code belows construct the Motzkin polynomial using [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl).

using Test #src
using DynamicPolynomials
@polyvar x y
motzkin = x^4*y^2 + x^2*y^4 + 1 - 3x^2*y^2

# The Motzkin polynomial is nonnegative but is not a sum of squares as we can verify numerically as follows.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.

using SumOfSquares
import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver)
@constraint(model, motzkin >= 0) # We constraint `motzkin` to be a sum of squares

optimize!(model)

# We see that the problem is detected as infeasible...

@test termination_status(model) == MOI.INFEASIBLE #src
termination_status(model)

# ... and that the dual solution is a certificate of the infeasibility of the problem.

@test dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE #src
dual_status(model)

# Even if the Motzkin polynomial is not a sum of squares, it can still be certified to be nonnegative using sums of squares.
# Indeed a polynomial is certified to be nonnegative if it is equal to a fraction of sums of squares.
# The Motzkin polynomial is equal to a fraction of sums of squares whose denominator is $x^2 + y^2$.
# This can be verified numerically as follows:

model = SOSModel(solver)
@constraint(model, (x^2 + y^2) * motzkin >= 0) # We constraint the `(x^2 + y^2) * motzkin` to be a sum of squares

optimize!(model)

# Now the problem is declared feasible by the solver...

@test termination_status(model) == MOI.OPTIMAL #src
termination_status(model)

# ... and the primal solution is a feasible point, hence it is a certificate of nonnegativity of the Motzkin polynomial.

@test primal_status(model) == MOI.FEASIBLE_POINT #src
primal_status(model)

# One may consider ourself lucky to have had the intuition that $x^2 + y^2$ would work as denominator.
# In fact, the search for the denominator can be carried out in parallel to the search of the numerator.
# In the example below, we search for a denominator with monomials of degrees from 0 to 2.
# If none is found, we can increase the maximum degree 2 to 4, 6, 8, ...
# This gives a hierarchy of programs to try in order to certify the nonnegativity of a polynomial by identifying it with a fraction of sum of squares polynomials.
# In the case of the Motzkin polynomial we now that degree 2 is enough since $x^2 + y^2$ works.

model = SOSModel(solver)
X = monomials([x, y], 0:2)

# We create a quadratic polynomial that is not necessarily a sum of squares since this is implied by the next constraint: `deno >= 1`.

@variable(model, deno, Poly(X))

# We want the denominator polynomial to be strictly positive, this prevents the trivial solution deno = 0 for instance.

@constraint(model, deno >= 1)
@constraint(model, deno * motzkin >= 0)
optimize!(model)

@test termination_status(model) == MOI.OPTIMAL #src
termination_status(model)

@test primal_status(model) == MOI.FEASIBLE_POINT #src
primal_status(model)

# We can check the denominator found by the program using `JuMP.value`

value(deno)

# Because a picture is worth a thousand words let's plot the beast.
# We can easily extend `Plots` by adding a recipe to plot bivariate polynomials.

using RecipesBase
@recipe function f(x::AbstractVector, y::AbstractVector, p::Polynomial)
    x, y, (x, y) -> p(variables(p) => [x, y])
end
import Plots
Plots.plot(
    range(-2, stop=2, length=100),
    range(-2, stop=2, length=100),
    motzkin,
    st = [:surface],
    seriescolor=:heat,
    colorbar=:none,
    clims = (-10, 80)
)
