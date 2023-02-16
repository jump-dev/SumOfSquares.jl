using DynamicPolynomials
@polyvar x[1:2]

a = 1
I = 0.05
ε = 0.08
γ = 1
F0 = [-x[2] - x[1] * (x[1] - 1) * (x[1] - a) + I, ε * (x[1] - γ * x[2])]

x1 = x2 = 0.0256
F = [f(x => [x[1] + x1, x[2] + x2]) for f in F0]

J = [j(x => zeros(2)) for j in differentiate(F, x)]

using LinearAlgebra
E = eigen(J)

λ = E.values[end]
w = E.vectors[:, end]

using SumOfSquares
r = 0.3
X = @set x[1]^2 + x[2]^2 ≤ r^2

import CSDP
model = SOSModel(CSDP.Optimizer)
@variable(model, γ)
@objective(model, Min, γ)
@variable(model, ϕN, Poly(monomials(x, 2:10)))
ϕ = w ⋅ x + ϕN
∇ϕ = differentiate(ϕ, x)
@constraint(model, -γ ≤ F ⋅ ∇ϕ - λ * ϕ, domain = X)
@constraint(model, F ⋅ ∇ϕ - λ * ϕ ≤ γ, domain = X)
optimize!(model)
solution_summary(model)

ϕ_opt = value(ϕ)

using Plots
x1s = x2s = range(-0.3, stop = 0.3, length = 40)
ϕs = ϕ_opt.(x1s', x2s)
contour(x1s, x2s, ϕs, levels=[-0.2, -0.1, 0, 0.1, 0.2, 0.5, 0.7])
θ = range(0, stop = 2π, length = 100)
plot!(r * cos.(θ), r * sin.(θ), label = "")
scatter!([0], [0], label = "")

using DifferentialEquations
function S(t, x1, x2, solver = DifferentialEquations.Tsit5())
    tspan = (0.0, t)
    prob = DifferentialEquations.ODEProblem((v, p, t) -> [f(x => v) for f in F], [x1, x2], tspan)
    traj = DifferentialEquations.solve(prob, solver, reltol=1e-4, abstol=1e-4)
    return traj[end]
end

using QuadGK
function laplace_average(f, λ, x1, x2, T = 10, args...)
    v, _ = quadgk(0, T, rtol=1e-3) do t
        s = S(t, x1, x2, args...)
        return f(S(t, x1, x2, args...)) * exp(-λ * t)
    end
    return v / T
end

lap(x1, x2) = laplace_average(v -> ϕ_opt(x => v), λ, x1, x2, 10)
laplace = lap.(x1s', x2s)

norm(laplace - ϕs)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

