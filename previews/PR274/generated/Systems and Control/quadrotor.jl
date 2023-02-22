using DynamicPolynomials
@polyvar x[1:6]
@polyvar u[1:2]
sinx5 = -0.166 * x[5]^3 + x[5]
cosx5 = -0.498 * x[5]^2 + 1
gn = 9.8
K = 0.89 / 1.4
d0 = 70
d1 = 17
n0 = 55
f = [
    x[3],
    x[4],
    0,
    -gn,
    x[6],
    -d0 * x[5] - d1 * x[6],
]
n_x = length(f)
g = [
    0         0
    0         0
    K * sinx5 0
    K * cosx5 0
    0         0
    0         n0
]
n_u = size(g, 2)

using SumOfSquares
rectangle = [1.7, 0.85, 0.8, 1, π/12, π/2, 1.5, π/12]
X = BasicSemialgebraicSet(FullSpace(), typeof(x[1] + 1.0)[])
for i in eachindex(x)
    addinequality!(X, x[i] + rectangle[i]) # x[i] >= -rectangle[i]
    addinequality!(X, rectangle[i] - x[i]) # x[i] <= rectangle[i]
end
X

using SparseArrays
x0 = zeros(n_x)
u0 = [gn/K, 0.0]

x_dot = f + g * u
A = map(differentiate(x_dot, x)) do a
    a(x => x0, u => u0)
end

B = map(differentiate(x_dot, u)) do b
    b(x => x0, u => u0)
end

import MatrixEquations
S, v, K = MatrixEquations.arec(A, B, 10, 100)

P, _, _ = MatrixEquations.arec(A - B * K, 0.0, 10.0)

nD = n_x + n_u
E = sparse(1:n_x, 1:n_x, ones(n_x), n_x, nD)
C = [A B]

using LinearAlgebra
import SCS
solver = optimizer_with_attributes(SCS.Optimizer, MOI.Silent() => true)
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

P = inv(Symmetric(value.(Q)))
using LinearAlgebra
F = cholesky(P)
K = -F.U[:, (n_x + 1):(nD)] \ F.U[:, 1:n_x] # That gives the following state feedback in polynomial form:

k = K * x

Px = inv(Symmetric(value.(Q[1:n_x, 1:n_x])))

Px = [I; K]' * P * [I; K]

eigen(Symmetric(Px * (A + B * K) + (A + B * K)' * Px)).values

function _create(model, d, P)
    if d isa Int
        return @variable(model, variable_type = P(monomials(x, 0:d)))
    else
        return d
    end
end
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

V

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
        @info("Solved in $(solve_time(model)) : γ ∈ ]$γ_min, $γ_max]")
    end
    if !isfinite(γ_max)
        error("Cannot find any infeasible γ")
    end
    return γ_min, k_best, s3_best
end

γ, k, s3 = γ_step(solver, V, γ, k, s3, [2, 2], 2, 1e-3, 10)

model, V = V_step(solver, V, γ, k, s3)
solution_summary(model)

V

γ, k, s3 = γ_step(solver, V, γ, k, s3, [2, 2], 2, 1e-3, 10)

model, V = V_step(solver, V, γ, k, s3)
solution_summary(model)

V

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

