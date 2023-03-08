# # Viability tube for quadrotor

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Systems and Control/quadrotor.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Systems and Control/quadrotor.ipynb)
# **Adapted from**: [YAP21, Section V.D] for the model defined in [B12], [M16, Section IV] and [M19, Section 6.1]
#
# [B12] Bouffard, Patrick.
# *On-board model predictive control of a quadrotor helicopter: Design, implementation, and experiments.*
# CALIFORNIA UNIV BERKELEY DEPT OF COMPUTER SCIENCES, 2012.
#
# [M16] Mitchell, Ian M., et al.
# *Ensuring safety for sampled data systems: An efficient algorithm for filtering potentially unsafe input signals.*
# 2016 IEEE 55th Conference on Decision and Control (CDC). IEEE, 2016.
#
# [M19] Mitchell, Ian M., Jacob Budzis, and Andriy Bolyachevets.
# *Invariant, viability and discriminating kernel under-approximation via zonotope scaling.*
# Proceedings of the 22nd ACM International Conference on Hybrid Systems: Computation and Control. 2019.
#
# [YAP21] Yin, H., Arcak, M., Packard, A., & Seiler, P. (2021).
# *Backward reachability for polynomial systems on a finite horizon.*
# IEEE Transactions on Automatic Control, 66(12), 6025-6032.

using Test #src
using DynamicPolynomials
@polyvar x[1:6]
@polyvar u[1:2]
sinx5 = -0.166 * x[5]^3 + x[5]
cosx5 = -0.498 * x[5]^2 + 1
gravity = 9.81
gain_u1 = 0.89 / 1.4
d0 = 70
d1 = 17
n0 = 55
f = [
    x[3],
    x[4],
    0,
    -gravity,
    x[6],
    -d0 * x[5] - d1 * x[6],
]
n_x = length(f)
g = [
    0               0
    0               0
    gain_u1 * sinx5 0
    gain_u1 * cosx5 0
    0               0
    0               n0
]
n_u = size(g, 2)

# The constraints below are the same as [YAP21, M16, M19] except
# [M16, M19] uses different bounds for `x[2]` and
# [M16] uses different bounds for `x[5]`

using SumOfSquares
rectangle = [1.7, 0.85, 0.8, 1, π/12, π/2, 1.5, π/12]
X = BasicSemialgebraicSet(FullSpace(), typeof(x[1] + 1.0)[])
for i in eachindex(x)
    lower = x[i] + rectangle[i] # x[i] >= -rectangle[i]
    upper = rectangle[i] - x[i] # x[i] <= rectangle[i]
    addinequality!(X, lower * upper) # -rectangle[i] <= x[i] <= rectangle[i]
end
X

## Controller for the linearized system

# The starting value for the Lyapunov function is the linear state-feedback
# that maintains the quadrotor at the origin [YAP21, Remark 3].

using SparseArrays
x0 = zeros(n_x)
u0 = [gravity/gain_u1, 0.0]

# The linearization of `f` is given by

x_dot = f + g * u
A = map(differentiate(x_dot, x)) do a
    a(x => x0, u => u0)
end

# The linearization of `g` is given by:

B = map(differentiate(x_dot, u)) do b
    b(x => x0, u => u0)
end

# We can compute the Linear-Quadratic Regulator using the same weight matrices
# as [YAP21](https://github.com/heyinUCB/Backward-Reachability-Analysis-and-Control-Synthesis)

import MatrixEquations
S, v, K = MatrixEquations.arec(A, B, 100, 10)
@test all(c -> real(c) < 0, v)

# The corresponding quadratic regulator is:

P, _, _ = MatrixEquations.arec(A - B * K, 0.0, 10.0)

# It corresponds to the following quadratic Lyapunov function:

V0 = x' * P * x

# ## γ-step

# It is a Lyapunov function for the linear system but not necessarily for the nonlinear system as well.
# We can however say that the γ-sublevel set `{x | x' P x ≤ γ}` is a (trivial) controlled invariant set for `γ = 0` (since it is empty).
# We can try to see if there a larger `γ` such that the γ-sublevel set is also controlled invariant using
# the γ step of [YAP21, Algorithm 1]

function _create(model, d, P)
    if d isa Int
        if P == SOSPoly   # `d` is the maxdegree while for `SOSPoly`,
            d = div(d, 2) # it is the monomial basis of the squares
        end
        return @variable(model, variable_type = P(monomials(x, 0:d)))
    else
        return d
    end
end

using LinearAlgebra
function base_model(solver, V, k, s3, γ)
    model = SOSModel(solver)
    V = _create(model, V, Poly)
    k = _create.(model, k, Poly)
    s3 = _create(model, s3, SOSPoly)
    ∂ = differentiate # Let's use a fancy shortcut
    @constraint(model, ∂(V, x) ⋅ (f + g * k) <= s3 * (V - γ)) # [YAP21, (E.2)]
    for r in inequalities(X) # `{V ≤ γ} ⊆ {r ≥ 0}` iff `r ≤ 0 => V ≥ γ`
        @constraint(model, V >= γ, domain = @set(r <= 0)) # [YAP21, (E.3)]
    end
    return model, V, k, s3
end

function γ_step(solver, V, γ_min, degree_k, degree_s3; γ_tol = 1e-1, max_iters = 10, γ_step = 0.5)
    γ0_min = γ_min
    γ_max = Inf
    num_iters = 0
    k_best = s3_best = nothing
    while true
        if γ_max - γ_min > γ_tol && num_iters < max_iters
            if isfinite(γ_max)
                γ = (γ_min + γ_max) / 2
            else
                γ = γ0_min + (γ_min - γ0_min + γ_step) * 2
            end
        elseif isnothing(k_best)
            @assert γ_min == γ0_min
            @info("Bisection finished without a feasible controller, let's find one")
            γ = γ0_min # Last run to compute a controller for the value of `γ` we know is feasible
        else
            break
        end
        model, V, k, s3 = base_model(solver, V, degree_k, degree_s3, γ)
        num_iters += 1
        @info("Iteration $num_iters/$max_iters : Solving with $(solver_name(model)) for `γ = $γ`")
        optimize!(model)
        @info("After $(solve_time(model)) seconds, terminated with $(termination_status(model)) ($(raw_status(model)))")
        if primal_status(model) == MOI.FEASIBLE_POINT
            @info("Feasible solution found")
            γ_min = γ
            k_best = value.(k)
            s3_best = value(s3)
        elseif dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
            @info("Infeasibility certificate found")
            if γ == γ0_min # This corresponds to the case above where we reached the tol or max iteration and we just did a last run at the value of `γ_min` provided by the user
                error("The value `$γ0_min` of `γ_min` provided is not feasible")
            end
            γ_max = γ
        else
            @warn("Giving up $(raw_status(model)), $(termination_status(model)), $(primal_status(model)), $(dual_status(model))")
            break
        end
        if γ != γ0_min
            @info("Refined interval : `γ ∈ [$γ_min, $γ_max[`")
        end
    end
    if !isfinite(γ_max) && max_iters > 0
        @warn("Cannot find any infeasible `γ` after $num_iters iterations")
    end
    return γ_min, k_best, s3_best
end

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
γ1, κ1, s3_1 = γ_step(solver, V0, 0.0, [2, 2], 2)
@test γ1 ≈ 0.5 atol = 1.e-1 #src

# Let's visualize now the controlled invariant set we have found:

using ImplicitPlots
using Plots
import ColorSchemes
function plot_lyapunovs(Vs, J; resolution = 1000, scheme = ColorSchemes.rainbow)
    xmax = rectangle[J[1]]
    ymax = rectangle[J[2]]
    rect_xs = [xmax, -xmax, -xmax, xmax, xmax]
    rect_ys = [ymax, ymax, -ymax, -ymax, ymax]
    p = plot(rect_xs, rect_ys, label="", color=:black)
    eliminated = x[setdiff(1:n_x, J)]
    for i in eachindex(Vs)
        V = subs(Vs[i], eliminated => zeros(length(eliminated)))
        linecolor = get(scheme, (i - 1) / max(length(Vs) - 1, 1))
        implicit_plot!(p, V; resolution, label="V$i", linecolor)
    end
    xlim = 1.05 * xmax
    ylim = 1.05 * ymax
    plot!(p, xlims=(-xlim, xlim), ylims=(-ylim, ylim), aspect_ratio = xmax / ymax)
    return p
end
Vs = [V0 - γ1]
plot_lyapunovs(Vs, [1, 2])

# ## V-step

# Let's now fix the control law that we have found and try to find a superset
# of the current controlled invariant set
# That corresponds to the V-step of [YAP21, Algorithm 1]:

_degree(d::Int) = d
_degree(V) = maxdegree(V)

function V_step(solver, V0, γ, k, s3)
    model, V, k, s3 = base_model(solver, _degree(V0), k, s3, γ)
    if !(V0 isa Int) # {V0 ≤ γ} ⊆ {V ≤ γ} iff V0 ≤ γ => V ≤ γ
        @constraint(model, V <= γ, domain = @set(V0 <= γ)) # [YAP21, (E.6)]
    end
    optimize!(model)
    return model, value(V)
end

model, V1 = V_step(solver, V0, γ1, κ1, s3_1)
solution_summary(model)

# We can see that the solver found a feasible solution.
# Let's compare it:

push!(Vs, V1 - γ1)
plot_lyapunovs(Vs, [1, 2])

# ## Continue iterating

# We could now find a better state feedback controller with this new Lyapunov.
# Given the plot, we see that `γ` cannot be much improved, we mostly
# want to find a better `κ` so let's set `max_iters` to zero

γ2, κ2, s3_2 = γ_step(solver, V0, γ1, [2, 2], 2, max_iters = 0)

# Let's see if this new controller allows to find a better Lyapunov.

model, V2 = V_step(solver, V0, γ2, κ2, s3_2)
solution_summary(model)

# It does not seem that we gained any improvement so let's stop:

push!(Vs, V2 - γ2)
plot_lyapunovs(Vs, [1, 2])
