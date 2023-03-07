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
@polyvar t
sinx5 = -0.166 * x[5]^3 + x[5]
cosx5 = -0.498 * x[5]^2 + 1
gravity = 9.81
K_const = 0.89 / 1.4
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
    K_const * sinx5 0
    K_const * cosx5 0
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
    addinequality!(X, x[i] + rectangle[i]) # x[i] >= -rectangle[i]
    addinequality!(X, rectangle[i] - x[i]) # x[i] <= rectangle[i]
end
X


# The starting value for `k` is the following linear state-feedback
# that maintains the quadrotor at the origin [YAP21, Remark 3].

using SparseArrays
x0 = zeros(n_x)
u0 = [gn/K, 0.0]

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

# The corresponding quadratic regulator is:

P, _, _ = MatrixEquations.arec(A - B * K, 0.0, 10.0)

V = x' * P * x

function _create(model, d, P)
    if d isa Int
        return @variable(model, variable_type = P(monomials(x, 0:d)))
    else
        return d
    end
end

using LinearAlgebra
function base_model(solver, V, k, s3, γ)
    T = 2;
    w = t*(T-t);
    model = SOSModel(solver)
    xt = [x; t]
    sos() = @variable(model, variable_type = SOSPoly(monomials(xt, 0:1)))
    function soseps()
        s = @variable(model, variable_type = Poly(monomials(xt, 0:2)))
        @constraint(model, s - 1e-4 in SOSCone())
        return s
    end
    s1 = @variable(model, variable_type = Poly(monomials(x, 0:2)))
    s2 = sos()
    s3 = sos()
    s4a = soseps()
    s4b = soseps()
    s4c = soseps()
    s4d = soseps()
    s4e = soseps()
    @variable(model, s4f, Poly(monomials(xt, 0:2))) # ?????
    s5a = sos()
    s5b = sos()
    s6a = sos()
    s6b = sos()
    s7a = sos()
    s7b = sos()
    s8a = sos()
    s8b = sos()
    s9a = sos()
    s9b = sos()
    s9c = sos()
    s9d = sos()
    s9e = sos()
    s9f = sos()
    @variable(model, k[1:2], Poly(monomials(xt, 0:2)))

    # dV/dt <= 0
    ∂ = differentiate # Let's use a fancy shortcut
    @constraint(model, ∂(V, x) ⋅ (f + g * k) + s2 * w <= s3 * (V - γ)) # [YAP21, (E.2)]
    # V(t,x)<=gamma implies rt<=0
    rt1 = (1.7 - x[1])*(x[1] + 1.7);
    rt2 = (0.85 - x[2])*(x[2] + 0.85);
    rt3 = (0.8 - x[3])*(x[3] + 0.8);
    rt4 = (1 - x[4])*(x[4] + 1);
    rt5 = (pi/12 - x[5])*(x[5] + pi/12);
    rt6 = (pi/2 - x[6])*(x[6] + pi/2);
    @constraint(model, s4a*rt1 + (V - γ) - s9a*w in SOSCone())
    @constraint(model, s4b*rt2 + (V - γ) - s9b*w in SOSCone())
    @constraint(model, s4c*rt3 + (V - γ) - s9c*w in SOSCone())
    @constraint(model, s4d*rt4 + (V - γ) - s9d*w in SOSCone())
    @constraint(model, s4e*rt5 + (V - γ) - s9e*w in SOSCone())
    @constraint(model, s4f*rt6 + (V - γ) - s9f*w in SOSCone())
    # V(t,x) <= gamma implies u <= uM
    uM = [1.5 + gravity/K_const; π/12];
    um = [-1.5 + gravity/K_const; -π/12];
    @constraint(model, uM[1] - k[1] + s5a*(V - γ) - s6a*w in SOSCone())
    @constraint(model, uM[2] - k[2] + s5b*(V - γ) - s6b*w in SOSCone())
    # V(t,x) <= gamma implies u >= um
    @constraint(model, k[1] - um[1] + s7a*(V - γ) - s8a*w in SOSCone())
    @constraint(model, k[2] - um[2] + s7b*(V - γ) - s8b*w in SOSCone())
    return model, V, k, s3
end

using MutableArithmetics
function γ_step(solver, V, γ_min, k_best, s3_best, degree_k, degree_s3, γ_tol, max_iters)
    γ0_min = γ_min
    γ_max = Inf
    num_iters = 0
    while γ_max - γ_min > γ_tol && num_iters < max_iters
        if isfinite(γ_max)
            γ = (γ_min + γ_max) / 2
        else
            γ = γ0_min + (γ_min - γ0_min + 1) * 2
        end
        model, V, k, s3 = base_model(solver, V, degree_k, degree_s3, γ)
        num_iters += 1
        @info("Iteration $num_iters/$max_iters : solving...")
        optimize!(model)
        if primal_status(model) == MOI.FEASIBLE_POINT
            γ_min = γ
            k_best = value.(k)
            s3_best = value(s3)
        elseif dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
            γ_max = γ
        else
            @warn("Giving up $(raw_status(model)), $(termination_status(model)), $(primal_status(model)), $(dual_status(model))")
            break
        end
        @info("Solved in $(solve_time(model)) : γ ∈ [$γ_min, $γ_max[")
    end
    if !isfinite(γ_max)
        error("Cannot find any infeasible γ")
    end
    return γ_min, k_best, s3_best
end


using MosekTools
solver = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
γ = 0.0
k = nothing
s3 = nothing
γ, k, s3 = γ_step(solver, V, γ, k, s3, [2, 2], 2, 1e-3, 10)

# This does not however take the constraints `X` into account.
# To take the constraint into account,
# we compute an ellipsoidal control invariant set using [LJ21, Corollary 9]
# For this, we first compute the descriptor system described in [LJ21, Proposition 5].
#
# [LJ21] Legat, Benoît, and Jungers, Raphaël M.
# *Geometric control of algebraic systems.*
# IFAC-PapersOnLine 54.5 (2021): 79-84.


nD = n_x + n_u
E = sparse(1:n_x, 1:n_x, ones(n_x), n_x, nD)
C = [A B]

# We know solve [LJ21, (13)]

model = Model(solver)
@variable(model, Q[1:nD, 1:nD] in PSDCone())
cref = @constraint(model, Symmetric(-C * Q * E' - E * Q * C') in PSDCone())
@constraint(model, rect_ref[i in 1:nD], Q[i, i] <= rectangle[i])
@variable(model, volume)
q = [Q[i, j] for j in 1:nD for i in 1:j]
@constraint(model, [volume; q] in MOI.RootDetConeTriangle(nD))
@objective(model, Max, volume)
optimize!(model)
solution_summary(model)

# We now have the control-Lyapunov function `V(x, u) = [x, u]' inv(Q) [x, u]`.
# In other words, The 1-sublevel set of `V(x, u)` is an invariant subset of `rectangle`
# with any state-feedback `κ(x)` such that `V(x, κ(x)) ≤ 1` for any `x` such that
# `min_u V(x, u) ≤ 1`.
# Such candidate `κ(x)` can therefore be chosen as `argmin_u V(x, u)`.
# Let `inv(Q) = U' * U` where `U = [Ux Uu]`. We have `V(x, u) = ||Ux * x + Uu * u||_2`.
# `κ(x)` is therefore the least-square solution of `Uu * κ(x) = -Ux * x`.
# This we find the linear state-feedback `κ(x) = K * x` where `K = -Uu \ Ux`.

P = inv(Symmetric(value.(Q)))
using LinearAlgebra
F = cholesky(P)
K = -F.U[:, (n_x + 1):(nD)] \ F.U[:, 1:n_x] # That gives the following state feedback in polynomial form:

# The corresponding polynomial form is given by:

k = K * x

# We now have two equivalent ways to obtain the Lyapunov function.
# Because `{V(x) ≤ 1} = {min_u V(x, u) ≤ 1}`,
# see the left-hand side as the projection of the ellipsoid on `x, u`.
# As the projection on the polar becomes simply cutting with the hyperplane `u = 0`,
# the polar of the projection is simply `Q[1:6, 1:6]` ! So

Px = inv(Symmetric(value.(Q[1:n_x, 1:n_x])))

# An alternative way is to use our linear state feedback.
# We know that `min_u V(x, u) = V(x, Kx)` so
Px = [I; K]' * P * [I; K]

# We can double check that this matrix is negative definite:

eigen(Symmetric(Px * (A + B * K) + (A + B * K)' * Px)).values

# Let's now find a valid Lyapunov function for the nonlinear system
# using that linear state feedback.
# That corresponds to the V-step of [YAP21, Algorithm 1]:

function base_model(solver, V, k, s3, γ)
    model = SOSModel(solver)
    V = _create(model, V, Poly)
    k = _create.(model, k, Poly)
    s3 = _create(model, s3, SOSPoly)
    ∂ = differentiate # Let's use a fancy shortcut
    @constraint(model, ∂(V, x) ⋅ (f + g * k) <= s3 * (V - γ)) # [YAP21, (E.2)]
    for r in inequalities(X)
        @constraint(model, V >= γ, domain = @set(r >= 0)) # [YAP21, (E.3)]
    end
    return model, V, k, s3
end

_degree(d::Int) = d
_degree(V) = maxdegree(V)

function V_step(solver, V0, γ, k, s3)
    model, V, k, s3 = base_model(solver, _degree(V0), k, s3, γ)
    if !(V0 isa Int)
        @constraint(model, V >= γ, domain = @set(V0 >= γ)) # [YAP21, (E.6)]
    end
    optimize!(model)
    return model, value(V)
end

γ = 1.0
s3 = 1.0
model, V = V_step(solver, 2, γ, k, s3)
solution_summary(model)

# The Lyapunov obtained is as follows

V

# We now try to find a state feedback that would improve γ


γ = 0.0
k = nothing
s3 = nothing
γ, k, s3 = γ_step(solver, V, γ, k, s3, [2, 2], 2, 1e-3, 10)

# We now try to find a new Lyapunov V:

model, V = V_step(solver, V, γ, k, s3)
solution_summary(model)

# The Lyapunov obtained is as follows

V

# We now try to improve γ again

γ, k, s3 = γ_step(solver, V, γ, k, s3, [2, 2], 2, 1e-3, 10)

# We now try to find a new Lyapunov V:

model, V = V_step(solver, V, γ, k, s3)
solution_summary(model)

# The Lyapunov obtained is as follows

V
