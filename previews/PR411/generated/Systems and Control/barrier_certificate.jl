using DynamicPolynomials
@polyvar x[1:2]

using SumOfSquares
using CSDP

solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

f = [ x[2],
     -x[1] + (1/3)*x[1]^3 - x[2]]

g₁ = -(x[1]+1)^2 - (x[2]+1)^2 + 0.16  # 𝒳ᵤ = {x ∈ R²: g₁(x) ≥ 0}

h₁ = -(x[1]-1.5)^2 - x[2]^2 + 0.25    # 𝒳₀ = {x ∈ R²: h₁(x) ≥ 0}

monos = monomials(x, 0:4)
@variable(model, B, Poly(monos))

ε = 0.001
@constraint(model, B >= ε, domain = @set(g₁ >= 0))

@constraint(model, B <= 0, domain = @set(h₁ >= 0))

using LinearAlgebra # Needed for `dot`
dBdt = dot(differentiate(B, x), f)
@constraint(model, -dBdt >= 0)

JuMP.optimize!(model)

JuMP.primal_status(model)

import OrdinaryDiffEq, Plots, ImplicitPlots
function phase_plot(f, B, g₁, h₁, quiver_scaling, Δt, X0, solver = OrdinaryDiffEq.Tsit5())
    X₀plot = ImplicitPlots.implicit_plot(h₁; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="X₀", linecolor=:blue)
    Xᵤplot = ImplicitPlots.implicit_plot!(g₁; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="Xᵤ", linecolor=:teal)
    Bplot  = ImplicitPlots.implicit_plot!(B; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="B = 0", linecolor=:red)
    Plots.plot(X₀plot)
    Plots.plot!(Xᵤplot)
    Plots.plot!(Bplot)
    ∇(vx, vy) = [fi(x[1] => vx, x[2] => vy) for fi in f]
    ∇pt(v, p, t) = ∇(v[1], v[2])
    function traj(v0)
        tspan = (0.0, Δt)
        prob = OrdinaryDiffEq.ODEProblem(∇pt, v0, tspan)
        return OrdinaryDiffEq.solve(prob, solver, reltol=1e-8, abstol=1e-8)
    end
    ticks = -5:0.5:5
    X = repeat(ticks, 1, length(ticks))
    Y = X'
    Plots.quiver!(X, Y, quiver = (x, y) -> ∇(x, y) / quiver_scaling, linewidth=0.5)
    for x0 in X0
        Plots.plot!(traj(x0), vars=(1, 2), label = nothing)
    end
    Plots.plot!(xlims = (-2, 3), ylims = (-2.5, 2.5))
end

phase_plot(f, value(B), g₁, h₁, 10, 30.0, [[x1, x2] for x1 in 1.2:0.2:1.7, x2 in -0.35:0.1:0.35])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
