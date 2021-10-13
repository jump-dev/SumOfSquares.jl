# # Polynomial Optimization

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/polynomial_optimization.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/polynomial_optimization.ipynb)
# **Contributed by**: Benoît Legat

# ## Introduction

# Consider the polynomial optimization problem [L09, Example 2.2] of
# minimizing the polynomial $x^3 - x^2 + 2xy - y^2 + y^3$
# over the polyhedron defined by the inequalities $x \ge 0, y \ge 0$ and $x + y \geq 1$.

# [L09] Lasserre, J. B.
# *Moments, positive polynomials and their applications*.
# World Scientific, **2009**.

using Test #src
using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
using SumOfSquares
S = @set x >= 0 && y >= 0 && x + y >= 1
p(x=>1, y=>0), p(x=>1//2, y=>1//2), p(x=>0, y=>1)

# The optimal solutions are $(x, y) = (1, 0)$ and $(x, y) = (0, 1)$ with objective value $0$ but [Ipopt](https://github.com/jump-dev/Ipopt.jl/) only finds the local minimum $(1/2, 1/2)$ with objective value $1/4$.

import Ipopt
model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
optimize!(model)
@test termination_status(model) == MOI.LOCALLY_SOLVED #src
@show termination_status(model)
@test value(a) ≈ 0.5 rtol=1e-5 #src
@show value(a)
@test value(b) ≈ 0.5 rtol=1e-5 #src
@show value(b)
@test objective_value(model) ≈ 0.25 rtol=1e-5 #src
@show objective_value(model)

# Note that the problem can be written equivalently as follows using [registered functions](https://jump.dev/JuMP.jl/stable/manual/nlp/#Register-a-function).

using Ipopt
model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
peval(a, b) = p(x=>a, y=>b)
register(model, :peval, 2, peval, autodiff=true)
@NLobjective(model, Min, peval(a, b))
optimize!(model)
@test termination_status(model) == MOI.LOCALLY_SOLVED #src
@show termination_status(model)
@test value(a) ≈ 0.5 rtol=1e-5 #src
@show value(a)
@test value(b) ≈ 0.5 rtol=1e-5 #src
@show value(b)
@test objective_value(model) ≈ 0.25 rtol=1e-5 #src
@show objective_value(model)

## Sum-of-Squares approach

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v0.21.6/installation/#Supported-solvers) for a list of the available choices.

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

# A Sum-of-Squares certificate that $p \ge \alpha$ over the domain `S`, ensures that $\alpha$ is a lower bound to the polynomial optimization problem.
# The following program searches for the largest lower bound and finds zero.

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c3, p >= α, domain = S)
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL #src
@show termination_status(model)
@test objective_value(model) ≈ 0.0 atol=1e-5 #src
@show objective_value(model)

# Using the solution $(1/2, 1/2)$ found by Ipopt of objective value $1/4$
# and this certificate of lower bound $0$ we know that the optimal objective value is in the interval $[0, 1/4]$
# but we still do not know what it is (if we consider that we did not try the solutions $(1, 0)$ and $(0, 1)$ as done in the introduction).
# If the dual of the constraint `c3` was atomic, its atoms would have given optimal solutions of objective value $0$ but that is not the case.

ν3 = moment_matrix(c3)
@test extractatoms(ν3, 1e-3) === nothing #src
extractatoms(ν3, 1e-3) # Returns nothing as the dual is not atomic

# Fortunately, there is a hierarchy of programs with increasingly better bounds that can be solved until we get one with atom dual variables.
# This comes from the way the Sum-of-Squares constraint with domain `S` is formulated.
# The polynomial $p - \alpha$ is guaranteed to be nonnegative over the domain `S` if there exists Sum-of-Squares polynomials $s_0$, $s_1$, $s_2$, $s_3$ such that
# ```math
# p - \alpha = s_0 + s_1 x + s_2 y + s_3 (x + y - 1).
# ```
# Indeed, in the domain `S`, $x$, $y$ and $x + y - 1$ are nonnegative so the right-hand side is a sum of squares hence is nonnegative.
# Once the degrees of $s_1$, $s_2$ and $s_3$ have been decided, the degree needed for $s_0$ will be determined but we have a freedom in choosing the degrees of $s_1$, $s_2$ and $s_3$.
# By default, they are chosen so that the degrees of $s_1 x$, $s_2 y$ and $s_3 (x + y - 1)$ match those of $p - \alpha$ but this can be overwritten using the `maxdegree` keyword argument.

# ### The maxdegree keyword argument

# The maximum total degree (i.e. maximum sum of the exponents of $x$ and $y$) of the monomials of $p$ is 3 so the constraint in the program above is equivalent to `@constraint(model, p >= α, domain = S, maxdegree = 3)`.
# That is, since $x$, $y$ and $x + y - 1$ have total degree 1, the sum of squares polynomials $s_1$, $s_2$ and $s_3$ have been chosen with maximum total degree $2$.
# Since these polynomials are sums of squares, their degree must be even so the next maximum total degree to try is 4.
# For this reason, the keywords `maxdegree = 4` and `maxdegree = 5` have the same effect in this example.
# In general, if the polynomials in the domain are not all odd or all even, each value of `maxdegree` has a different effect in the choice of the maximum total degree of some $s_i$.

function sos(deg)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = S, ideal_certificate = Certificate.MaxDegree(SOSCone(), MonomialBasis, 4), maxdegree = deg)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL #src
    @show termination_status(model)
    @test objective_value(model) ≈ 0.0 atol=1e-5 #src
    @show objective_value(model)
    display(lagrangian_multipliers(c))
    return c
end
model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c5, p >= α, domain = S, maxdegree = 5)
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL #src
@show termination_status(model)
@test objective_value(model) ≈ 0.0 atol=1e-5 #src
@show objective_value(model)

# As shown below, for `maxdegree = 5`, the dual variable is atomic as it is the moments of the measure
# $$0.5 \delta(x-1, y) + 0.5 \delta(x, y-1)$$
# where $\delta(x, y)$ is the dirac measure centered at $(0, 0)$.
# Therefore the program provides both a certificate that $0$ is a lower bound and a certificate that it is also an upper bound since it is attained at the global minimizers $(1, 0)$ and $(0, 1)$.

ν5 = moment_matrix(c5)
atoms5 = extractatoms(ν5, 1e-3) #src
@test atoms5.atoms[1].weight ≈ 0.5 rtol=1e-2 #src
@test atoms5.atoms[2].weight ≈ 0.5 rtol=1e-2 #src
@test atoms5.atoms[1].center[2:-1:1] ≈ atoms5.atoms[2].center[1:2] rtol=1e-2 #src
extractatoms(ν5, 1e-3)

# ## A deeper look into atom extraction

# The `extractatoms` function requires a `ranktol` argument that we have set to `1e-3` in the preceding section.
# This argument is used to transform the dual variable into a system of polynomials equations whose solutions give the atoms.
# This transformation uses the SVD decomposition of the moment matrix and discards the equations corresponding to a singular value lower than `ranktol`.
# When this system of equation has an infinite number of solutions, `extractatoms` concludes that the measure is not atomic.
# For instance, with `maxdegree = 3`, we obtain the system
# $$x + y = 1$$
# which contains a whole line of solution.
# This explains `extractatoms` returned `nothing`.

ν3 = moment_matrix(c3)
SumOfSquares.MultivariateMoments.computesupport!(ν3, 1e-3)
@test length(ν3.support.I.p) == 1 #src

# With `maxdegree = 5`, we obtain the system
# ```math
# \begin{aligned}
#   x + y & = 1\\
#   y^2 & = y\\
#   xy & = 0\\
#   x^2 + y & = 1
# \end{aligned}
# ```

ν5 = moment_matrix(c5)
SumOfSquares.MultivariateMoments.computesupport!(ν5, 1e-3)

# This system can be reduced to the equivalent system
# ```math
# \begin{aligned}
#   x + y & = 1\\
#   y^2 & = y
# \end{aligned}
# ```
# which has the solutions $(0, 1)$ and $(1, 0)$.

SemialgebraicSets.computegröbnerbasis!(ideal(ν5.support))
ν5.support
@test length(ν5.support.I.p) == 2 #src

# The function `extractatoms` then reuses the matrix of moments to find the weights $1/2$, $1/2$ corresponding to the diracs centered respectively at $(0, 1)$ and $(1, 0)$.
# This details how the function obtained the result
# $$0.5 \delta(x-1, y) + 0.5 \delta(x, y-1)$$
# given in the previous section.

# ## HomotopyContinuation

# As discussed in the previous section, the atom extraction relies on the solution
# of a system of algebraic equations. The `extractatoms` function takes an optional
# `algebraic_solver` argument that is used to solve this system of equation.
# If no solver is provided, the default solver of SemialgebraicSets.jl is used which
# currently computes the Gröbner basis, then the multiplication matrices and
# then the Schur decomposition of a random combination of these matrices.
# As the system of equations is obtained from a numerical solution and is represented
# using floating point coefficients, homotopy continuation is recommended as it is
# more numerically robust than Gröbner basis computation.
# The following uses homotopy continuation to solve the system of equations.

using HomotopyContinuation
algebraic_solver = SemialgebraicSetsHCSolver(; excess_residual_tol = 2e-2, real_tol = 2e-2, compile = false)
atoms5 = extractatoms(ν5, 1e-3, algebraic_solver) #src
@test length(atoms5.atoms) == 2 #src
@test atoms5.atoms[1].weight + atoms5.atoms[2].weight ≈ 1.0 rtol=1e-2 #src
@test atoms5.atoms[1].center[2:-1:1] ≈ atoms5.atoms[2].center[1:2] rtol=1e-2 #src
extractatoms(ν5, 1e-3, algebraic_solver)

# As the system has 3 equations for 2 variables and the coefficients of the equations
# are to be treated with tolerance since they originate from the solution of an SDP,
# we need to set `excess_residual_tol` and `real_tol` to a high tolerance otherwise,
# HomotopyContinuation would consider that there is no solution.
# Indeed, as the system is overdetermined (it has more equations than variables)
# HomotopyContinuation expects to have excess solution hence it filters out
# excess solution among the solution found. It determines which solution are in excess
# by comparing the infinity norm of the residuals of the equations at the solution with `excess_residual_tol`.
# It also filters out solution for which the absolute value of the imaginary part of one of the entry
# is larger than `real_tol` and strips out the imaginary part.
# The raw solutions obtained by HomotopyContinuation can be obtained as follows:

F = HomotopyContinuation.System(ν5.support)
res = HomotopyContinuation.solve(F, algebraic_solver.options...)
r = path_results(res) #src
@test length(r) == 4 #src
@test all(HomotopyContinuation.is_excess_solution, r) #src
path_results(res)

# The printed `residual` above shows why `2e-2` allows to filter how the 2 actual
# solutions from the 2 excess solutions.

# ## Truncated moment series
#
# A third alternative is implemented in https://github.com/bmourrain/MultivariateSeries.jl
# As it is not released yet, we need to add the package as follows:

using Pkg
pkg"add https://github.com/blegat/MultivariateSeries.jl#ds18"

# We will use `import`, not `using` to make it clear which function is in which
# package:

const MM = MultivariateMoments
import MultivariateSeries
const MS = MultivariateSeries

function decompose_truncation(ν::MM.MomentMatrix)
    μ = MM.measure(ν)
    s = MultivariateSeries.series(
        [monomial(moment) => moment_value(moment) for moment in MM.moments(μ)]
    )
    t = MultivariateSeries.truncate(s, maxdegree(s) - 1)
    weights, sols = MultivariateSeries.decompose(t)
    centers = [sols[:, i] for i in 1:size(sols, 2)]
    return AtomicMeasure(
        variables(ν),
        WeightedDiracMeasure.(centers, weights),
    )
end

# Let's try it first for `ν3`:

η3 = decompose_truncation(ν3)

# Note that the moments do not match (which is to be expected as the `(0.5, 0.5)` is not a minimizer).
# We can verify this by comparing the moments in the moment matrix:

μ3 = MM.measure(ν3)

# with the corresponding moments of the atomic measure:

MM.measure(η3, monomials(μ3))

# We try it out for `ν5` now:

η5 = decompose_truncation(ν5)

# Note the accuracy of the solution. This is thanks to the truncation that
# discards the highest degree moments which are less accurate.
# We can verify that the moments indeed match this time:

μ5 = MM.measure(ν5)

# with the corresponding moments of the atomic measure:

MM.measure(η5, monomials(μ5))

# Note that while the accuracy isn't so good for the first 5 moment values,
# it is quite accurate for the last ones.
# This is to be expected as these first 5 moments corresponds to the degree 4
# monomials which were truncated.
