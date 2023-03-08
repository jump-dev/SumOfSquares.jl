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

using SumOfSquares
rectangle = [1.7, 0.85, 0.8, 1, π/12, π/2, 1.5, π/12]
X = BasicSemialgebraicSet(FullSpace(), typeof(x[1] + 1.0)[])
for i in eachindex(x)
    lower = x[i] + rectangle[i] # x[i] >= -rectangle[i]
    upper = rectangle[i] - x[i] # x[i] <= rectangle[i]
    addinequality!(X, lower * upper) # -rectangle[i] <= x[i] <= rectangle[i]
end
X

# Controller for the linearized system

using SparseArrays
x0 = zeros(n_x)
u0 = [gravity/gain_u1, 0.0]

x_dot = f + g * u
A = map(differentiate(x_dot, x)) do a
    a(x => x0, u => u0)
end

B = map(differentiate(x_dot, u)) do b
    b(x => x0, u => u0)
end

import MatrixEquations
S, v, K = MatrixEquations.arec(A, B, 100, 10)

P, _, _ = MatrixEquations.arec(A - B * K, 0.0, 10.0)

V0 = x' * P * x

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

push!(Vs, V1 - γ1)
plot_lyapunovs(Vs, [1, 2])

γ2, κ2, s3_2 = γ_step(solver, V0, γ1, [2, 2], 2, max_iters = 0)

model, V2 = V_step(solver, V0, γ2, κ2, s3_2)
solution_summary(model)

push!(Vs, V2 - γ2)
plot_lyapunovs(Vs, [1, 2])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

