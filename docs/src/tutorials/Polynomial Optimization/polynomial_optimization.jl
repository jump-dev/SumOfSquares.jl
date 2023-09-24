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

# ## Local search

# A local solver only uses the **local** information given by the the value, gradient and hessian
# of the objective function and constraints at a given solution. When it converges, it therefore only
# guarantees that the found solution is a **local** minimum.
# In this example, the optimal solutions are $(x, y) = (1, 0)$ and $(x, y) = (0, 1)$ with objective value $0$ but
# [Ipopt](https://github.com/jump-dev/Ipopt.jl/) only finds the local minimum $(1/2, 1/2)$ with objective value $1/4$.

import Ipopt
model = Model(Ipopt.Optimizer)
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@objective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
optimize!(model)

# As we can see below, the termination status is `LOCALLY_SOLVED` and not of `OPTIMAL`
# because Ipopt only guarantees **local** optimality.

@test termination_status(model) == MOI.LOCALLY_SOLVED #src
@test objective_value(model) ≈ 0.25 rtol=1e-5 #src
solution_summary(model)

# Indeed, the solution found is not globally optimal:

@test value(a) ≈ 0.5 rtol=1e-5 #src
@test value(b) ≈ 0.5 rtol=1e-5 #src
value(a), value(b)

# Note that the problem can be written equivalently as follows using [registered functions](https://jump.dev/JuMP.jl/stable/manual/nlp/#Register-a-function).
# The difference is that the gradient and hessian will be computed via the *Symbolic Differentiation* provided
# by MultivariatePolynomials instead of JuMP's *Automatic Differentiation*:

f(a, b) = p(x => a, y => b)
∇p = differentiate(p, [x, y])
function ∇f(g, a, b)
    for i in eachindex(g)
        g[i] = ∇p[i](x => a, y => b)
    end
end
∇²p = differentiate(∇p, [x, y])
function ∇²f(H, a, b)
    for j in axes(∇²p, 2)
        for i in j:size(∇²p, 1)
            H[i, j] = ∇²p[i, j](x => a, y => b)
        end
    end
end
using Ipopt
gmodel = Model(Ipopt.Optimizer)
@variable(gmodel, a >= 0)
@variable(gmodel, b >= 0)
@constraint(gmodel, a + b >= 1)
register(gmodel, :f, 2, f, ∇f, ∇²f)
@NLobjective(gmodel, Min, f(a, b))
optimize!(gmodel)

# Even if we have the algebraic expressions of gradient and hessian,
# Ipopt is not using these symbolic expressions but only local information
# hence it can still only provide local guarantees:

@test termination_status(gmodel) == MOI.LOCALLY_SOLVED #src
@test objective_value(gmodel) ≈ 0.25 rtol=1e-5 #src
solution_summary(gmodel)

# and the same solution is found:

@test value(a) ≈ 0.5 rtol=1e-5 #src
@test value(b) ≈ 0.5 rtol=1e-5 #src
value(a), value(b)

# ## QCQP approach

import Alpine, HiGHS, Ipopt, Pavito
ipopt = optimizer_with_attributes(
    Ipopt.Optimizer,
    MOI.Silent() => true,
)
highs = optimizer_with_attributes(
    HiGHS.Optimizer,
    "presolve" => "on",
    "log_to_console" => false,
)
pavito = optimizer_with_attributes(
    Pavito.Optimizer,
    MOI.Silent() => true,
    "mip_solver" => highs,
    "cont_solver" => ipopt,
    "mip_solver_drives" => false,
)
alpine = optimizer_with_attributes(
    Alpine.Optimizer,
    "nlp_solver" => ipopt,
    "mip_solver" => pavito,
)
set_optimizer(model, () -> PolyJuMP.QCQP.Optimizer(MOI.instantiate(alpine)))
optimize!(model)

# We can see that it found the optimal solution

termination_status(model), value(a), value(b)

# ## Sum-of-Squares approach

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.

import SCS
scs = SCS.Optimizer
import Dualization
dual_scs = Dualization.dual_optimizer(scs)


# A Sum-of-Squares certificate that $p \ge \alpha$ over the domain `S`, ensures that $\alpha$ is a lower bound to the polynomial optimization problem.
# The following program searches for the largest lower bound and finds zero.

model = SOSModel(dual_scs)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c3, p >= α, domain = S)
optimize!(model)

# This time, the termination status is `OPTIMAL` but this does not necessarily mean that we found
# the optimal solution to the polynomial optimization problem.
# This only means that CSDP founds an optimal solution to the Sum-of-Squares relaxation.

@test termination_status(model) == MOI.OPTIMAL #src
@test objective_value(model) ≈ 0.0 atol=1e-3 #src
solution_summary(model)

# The feasibility of the primal solution guarantees that the objective value `0` is a lower bound
# to the polynomial optimization problem.
# The optimality means that it's the best lower bound we can get (at this degree of the hierarcy).
# Using the solution $(1/2, 1/2)$ found by Ipopt of objective value $1/4$
# and this certificate of lower bound $0$ we know that the optimal objective value is in the interval $[0, 1/4]$
# but we still do not know what it is (if we consider that we did not try the solutions $(1, 0)$ and $(0, 1)$ as done in the introduction).
# If the dual of the constraint `c3` was atomic, its atoms would have given optimal solutions of objective value $0$ but that is not the case.

ν3 = moment_matrix(c3)
@test atomic_measure(ν3, 1e-3) === nothing #src
atomic_measure(ν3, 1e-3) # Returns nothing as the dual is not atomic

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

model = SOSModel(dual_scs)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c4, p >= α, domain = S, maxdegree = 4)
optimize!(model)

# We can see that the basis of the moment matrix didn't increase:

@test length(moment_matrix(c4).basis.monomials) == 3 #src
moment_matrix(c4)

# This is because of the Newton polytope reduction that determined that gram matrix will
# be zero for these degrees so it reduced the problem back to the equivalent of `maxdegree` 3
# Let's turn this off with `newton_polytope = nothing`

function sos(solver, deg)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = S, maxdegree = deg, newton_polytope = nothing)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL #src
    @test objective_value(model) ≈ 0.0 atol=1e-3 #src
    return model
end
dual_model4 = sos(dual_scs, 4)
nothing #hide

# We see that the lower bound is still 0:

solution_summary(dual_model4)

# Let's now look at which solution we can extract from the moment matrix:

dual_ν4 = moment_matrix(dual_model4[:c])

# Looking at the singular values, `4` seems to be a reasonable rank:

using LinearAlgebra
svdvals(Matrix(dual_ν4.Q))

# The solution we extract is `(0.5, 0.5)` which is the solution found by Ipopt:

dual_atoms4 = atomic_measure(dual_ν4, FixedRank(4)) #src
@test dual_atoms4.atoms[1].center ≈ [0.5, 0.5] rtol=1e-1 #src
atomic_measure(dual_ν4, FixedRank(4))

# This process is quite sensitive numerically so let's try to solve it without dualization as well:

model4 = sos(scs, 4)
nothing #hide

# We see that the lower bound is again 0:

solution_summary(model4)

# The moment matrix is the following

ν4 = moment_matrix(model4[:c])

# Looking at the singular values, `3` seems to be a reasonable rank:

svdvals(Matrix(ν4.Q))

# This time, the dual variable is atomic as it is the moments of the measure
# $$0.5 \delta(x-1, y) + 0.5 \delta(x, y-1)$$
# where $\delta(x, y)$ is the dirac measure centered at $(0, 0)$.
# Therefore the program provides both a certificate that $0$ is a lower bound and a certificate that it is also an upper bound since it is attained at the global minimizers $(1, 0)$ and $(0, 1)$.

atoms4 = atomic_measure(ν4, FixedRank(3)) #src
@test atoms4.atoms[1].center[2:-1:1] ≈ atoms4.atoms[2].center[1:2] rtol=1e-1 #src
atomic_measure(ν4, FixedRank(3))

# ## A deeper look into atom extraction

# The moment matrix is transformed into a system of polynomials equations whose solutions give the atoms.
# This transformation uses the SVD decomposition of the moment matrix and discards the equations corresponding to the lowest singular values.
# When this system of equation has an infinite number of solutions, `atomic_measure` concludes that the measure is not atomic.
# For instance, with `maxdegree = 3`, we obtain the system
# $$x + y = 1$$
# which contains a whole line of solution.
# This explains `atomic_measure` returned `nothing`.

ν3 = moment_matrix(c3)
SumOfSquares.MultivariateMoments.compute_support!(ν3, LeadingRelativeRankTol(1e-3))
@test length(ν3.support.I.p) == 1 #src

# With `maxdegree = 4`, we obtain the system
# ```math
# \begin{aligned}
#   x + y & = 1\\
#   y^2 & = y\\
#   xy & = 0\\
#   -y + y^2 - x*y & = 0
#   y^2 - 2y + 1 & = x^2
# \end{aligned}
# ```

ν4 = moment_matrix(model4[:c])
SumOfSquares.MultivariateMoments.compute_support!(ν4, FixedRank(3))

# This system can be reduced to the equivalent system
# ```math
# \begin{aligned}
#   x + y & = 1\\
#   y^2 & = y
# \end{aligned}
# ```
# which has the solutions $(0, 1)$ and $(1, 0)$.

SemialgebraicSets.compute_gröbner_basis!(ideal(ν4.support))
ν4.support
@test length(ν4.support.I.p) == 2 #src

# The solutions of this system then give the minimizers

@test length(collect(ν4.support)) == 2 #src
collect(ν4.support)

# The function `atomic_measure` then reuses the matrix of moments to find the weights $1/2$, $1/2$ corresponding to the diracs centered respectively at $(0, 1)$ and $(1, 0)$.
# This details how the function obtained the result
# $$0.5 \delta(x-1, y) + 0.5 \delta(x, y-1)$$
# given in the previous section.

# ## HomotopyContinuation

# As discussed in the previous section, the atom extraction relies on the solution
# of a system of algebraic equations. The `atomic_measure` function takes an optional
# `algebraic_solver` argument that is used to solve this system of equation.
# If no solver is provided, the default solver of SemialgebraicSets.jl is used which
# currently computes the Gröbner basis, then the multiplication matrices and
# then the Schur decomposition of a random combination of these matrices.
# As the system of equations is obtained from a numerical solution and is represented
# using floating point coefficients, homotopy continuation is recommended as it is
# more numerically robust than Gröbner basis computation.
# The following uses homotopy continuation to solve the system of equations.

using HomotopyContinuation
algebraic_solver = SemialgebraicSetsHCSolver(; excess_residual_tol = 1e-1, real_tol = 1e-1, compile = false)
atoms4 = atomic_measure(ν4, FixedRank(3), Echelon(), algebraic_solver) #src
@test length(atoms4.atoms) == 2 #src
@test atoms4.atoms[1].weight + atoms4.atoms[2].weight ≈ 1.0 rtol=1e-1 #src
@test atoms4.atoms[1].center[2:-1:1] ≈ atoms4.atoms[2].center[1:2] rtol=1e-1 #src
atomic_measure(ν4, FixedRank(3), Echelon(), algebraic_solver)

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

F = HomotopyContinuation.System(ν4.support)
res = HomotopyContinuation.solve(F, algebraic_solver.options...)
r = path_results(res) #src
@test length(r) == 4 #src
@test all(HomotopyContinuation.is_excess_solution, r) #src
path_results(res)

# The printed `residual` above shows why `1e-1` allows to filter how the 2 actual
# solutions from the 2 excess solutions.
