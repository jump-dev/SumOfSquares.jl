using DynamicPolynomials
using LinearAlgebra
@polyvar x[1:2]

a = 1
I = 0.05
ε = 0.08
γ = 1
F0 = [-x[2] - x[1] * (x[1] - 1) * (x[1] - a) + I, ε * (x[1] - γ * x[2])]

x_eq = let z = 0.0256
    for _ in 1:5
        z -= (-z - z * (z - 1)^2 + I) / (-2 + 4z - 3z^2)
    end
    z
end
F = [f(x => x .+ x_eq) for f in F0]

J = [j(x => zeros(2)) for j in differentiate(F, x)]

E = eigen(J')
λ = E.values[end]
w = E.vectors[:, end]

if w[2] < 0
    w = -w
end

using SumOfSquares
r = 0.3
@polyvar y[1:2]
G = [f(x => r .* y) / r for f in F]
X = @set y[1]^2 + y[2]^2 ≤ 1

function fitzhugh_nagumo(solver, N, M)
    model = SOSModel(solver)
    @variable(model, η)
    @objective(model, Min, η)
    @variable(model, ϕN, Poly(monomials(y, 2:N)))
    ϕ = r * (w ⋅ y) + ϕN
    residual = G ⋅ differentiate(ϕ, y) - λ * ϕ
    @constraint(model, -η ≤ residual, domain = X, maxdegree = M + 2)
    @constraint(model, residual ≤ η, domain = X, maxdegree = M + 2)
    optimize!(model)
    return ϕ(y => x ./ r), model
end

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
ϕ, model = fitzhugh_nagumo(solver, 10, 20)
solution_summary(model)

ϕ_opt = value(ϕ)
residual_opt = F ⋅ differentiate(ϕ_opt, x) - λ * ϕ_opt

using Plots
x1s = x2s = range(-r, r, length = 40)
inside = [u^2 + v^2 ≤ r^2 for v in x2s, u in x1s]
ϕs = [
    inside[i, j] ? ϕ_opt(x => [u, v]) : NaN for
    (i, v) in enumerate(x2s), (j, u) in enumerate(x1s)
]
contour(
    x1s,
    x2s,
    ϕs;
    levels = [-0.2, -0.1, 0, 0.1, 0.2, 0.5, 0.7],
    aspect_ratio = :equal,
    xlabel = "x₁",
    ylabel = "x₂",
)
θ = range(0, 2π, length = 100)
plot!(r * cos.(θ), r * sin.(θ), label = "")
scatter!([0], [0], label = "")

import OrdinaryDiffEqTsit5
using QuadGK
function laplace_average(f, λ, x1, x2, T = 10.0)
    prob = OrdinaryDiffEqTsit5.ODEProblem(
        (v, p, t) -> [g(x => v) for g in F],
        [x1, x2],
        (0.0, T),
    )
    traj = OrdinaryDiffEqTsit5.solve(
        prob,
        OrdinaryDiffEqTsit5.Tsit5();
        reltol = 1e-11,
        abstol = 1e-13,
    )
    v, _ = quadgk(t -> f(traj(t)) * exp(-λ * t), 0.0, T; rtol = 1e-7)
    return v / T
end

lap(x1, x2) = laplace_average(v -> ϕ_opt(x => v), λ, x1, x2)
laplace = [
    inside[i, j] ? lap(u, v) : NaN for
    (i, v) in enumerate(x2s), (j, u) in enumerate(x1s)
]

maximum(abs, laplace[inside] - ϕs[inside])

inner = [u^2 + v^2 ≤ (r / 2)^2 for v in x2s, u in x1s]
maximum(abs, laplace[inner] - ϕs[inner])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
