# # Getting started

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started.ipynb)
# **Adapted from**: SOSTOOLS' SOSDEMO1 (See Section 4.1 of [SOSTOOLS User's Manual](http://sysos.eng.ox.ac.uk/sostools/sostools.pdf)) and Example 2.4 of [PJ08]
#
# P. Parrilo and A. Jadbabaie
# *Approximation of the joint spectral radius using sum of squares*.
# Linear Algebra and its Applications, Elsevier (2008), 428, 2385-2402

using Test #src
using DynamicPolynomials
@polyvar x y
p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4

# We need to pick an SDP solver, see [here](http://jump.dev/JuMP.jl/dev/installation/#Getting-Solvers-1) for a list of the available choices.
# We use `SOSModel` instead of `Model` to be able to use the `>=` syntax for Sum-of-Squares constraints.

using SumOfSquares
import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver)
con_ref = @constraint(model, p >= 0)
optimize!(model)
@test primal_status(model) == MOI.FEASIBLE_POINT #src
primal_status(model)

# We see above that the solver found a feasible solution.
# We now inspect this solution:

q = gram_matrix(con_ref)
Q = getmat(q) #src
@test isapprox(Q[1, 1], 2, rtol=1e-5) #src
@test isapprox(Q[1, 2], 1, rtol=1e-5) #src
@test isapprox(Q[3, 3], 5, rtol=1e-5) #src
@test abs(Q[2, 3]) < 1e-5 #src
@test isapprox(Q[2, 2] + 2Q[1, 3], -1, rtol=1e-5) #src

# We can get the SOS decomposition from the gram matrix as follows:

sosdec = SOSDecomposition(q)
@test isapprox(sosdec, sos_decomposition(con_ref)) #src
@test isapprox(sum(sosdec.ps.^2), p; rtol=1e-4, ztol=1e-6) #src

# We now seek for the SOS decomposition of the following polynomial:

p = 4*x^4*y^6 + x^2 - x*y^2 + y^2

# We build the same model as previously with this new polynomial.
# Here we can use `Model` instead of `SOSModel` as we explicitly constrain
# `p` to belong to the SOS cone with `p in SOSCone()`.

model = Model(solver)
con_ref = @constraint(model, p in SOSCone())
optimize!(model)
@test primal_status(model) == MOI.FEASIBLE_POINT #src
primal_status(model)

# We can query the SOS decomposition direction from the constraint reference
# as follows:

sos_decomposition(con_ref)
# p should be GramMatrix([1, 0, -1/2, 0, -1, 1, 0, -2/3, 0, 4/3, 0, 0, 2, 0, 4], [y, x, x*y, x*y^2, x^2*y^3]) #src
