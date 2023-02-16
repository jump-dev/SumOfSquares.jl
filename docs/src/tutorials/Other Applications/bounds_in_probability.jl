# # Bounds in Probability

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Other Applications/bounds_in_probability.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Other Applications/bounds_in_probability.ipynb)
# **Adapted from**: SOSTOOLS' SOSDEMO8 (See Section 4.8 of [SOSTOOLS User's Manual](http://sysos.eng.ox.ac.uk/sostools/sostools.pdf))

using Test #src

# The probability adds up to one.

μ0 = 1

# The mean is one.

μ1  = 1

# The standard deviation is 1/2.

σ = 1/2

# The second moment `E(x^2)` is:

μ2 = σ^2 + μ1^2

# We define the moments as follows:

using DynamicPolynomials
@polyvar x
monos = [1, x, x^2]
using SumOfSquares
μ = measure([μ0, μ1, μ2], monos)

# We need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.8/installation/#Supported-solvers) for a list of the available choices.
# We use `SOSModel` instead of `Model` to be able to use the `>=` syntax for Sum-of-Squares constraints.

using CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

# We create a polynomial with the monomials in `monos` and
# JuMP decision variables as coefficients as follows:

@variable(model, poly, Poly(monos))

# Nonnegative on the support:

K = @set 0 <= x && x <= 5
con_ref = @constraint(model, poly >= 0, domain = K)

# Greater than one on the event:

@constraint(model, poly >= 1, domain = (@set 4 <= x && x <= 5))

# The bound (we use `LinearAlgebra` for the `⋅` syntax for the scalar product):

using LinearAlgebra
@objective(model, Min, poly ⋅ μ)

# We verify that we found a feasible solution:

optimize!(model)
@test primal_status(model) == MOI.FEASIBLE_POINT #src
primal_status(model)

# The objective value is `1/37`:

@test isapprox(objective_value(model), 1/37, rtol=1e-5) #src
objective_value(model)

# The solution is `(12x-11)^2 / 37^2`:

@test sos_decomposition(con_ref, K, 1e-4) isa SOSDecompositionWithDomain #src
@test isapprox(value(poly), ((12/37)x-11/37)^2, rtol=1e-3) #src
value(poly) * 37^2
