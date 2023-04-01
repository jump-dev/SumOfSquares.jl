# # Barrier certificates for collision avoidance

# **Adapted from**: Example 2 of "Engineering applications of SOS polynomials" by Georgina Hall, from 
# the 2019 AMS Short Course on Sum of Squares: Theory and Applications

using Test #src
using DynamicPolynomials
@polyvar x[1:2]

using SumOfSquares
using CSDP

# We need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.8/installation/#Supported-solvers) for a list of the available choices.
# We use `SOSModel` instead of `Model` to be able to use the `>=` syntax for Sum-of-Squares constraints.
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

# We define below the vector field ``\text{d}x/\text{d}t = f``

f = [ x[2],
     -x[1] + (1/3)*x[1]^3 - x[2]]

# We now define the functions that describes the initial and unsafe sets

g₁ = -(x[1]+1)^2 - (x[2]+1)^2 + 0.16
h₁ = -(x[1]-1.5)^2 - x[2]^2 + 0.25

# Define SOS barrier function B
monos = monomials(x, 0:4)
@variable(model, B, Poly(monos))

# Define barrier certificate constraints 

# B(x) = ε + σ₀(x) + σ₁(x)g₁(x) where σᵢ are SOS and ε > 0
ε = 0.001
@variable(model, σ₁, Poly(monos))
@constraint(model, B - ε - σ₁*g₁ >= 0)

# -B(x) = τ₀(x) + τ₁(x)h₁(x) where τᵢ are SOS
@variable(model, τ₁, Poly(monos))
@constraint(model, -B - τ₁*h₁ >= 0)

# -Ḃ is SOS
using LinearAlgebra # Needed for `dot`
dBdt = dot(differentiate(B, x), f)
@constraint(model, -dBdt >= 0)

# The model is ready to be optimized by the solver:
JuMP.optimize!(model)

# We verify that the solver has found a feasible solution:
JuMP.primal_status(model)
@test JuMP.primal_status(model) == MOI.FEASIBLE_POINT #src

# Plot the phase plot with the 0-level set of the barrier function, and the boundary of the initial and unsafe sets
import DifferentialEquations, Plots, ImplicitPlots
function phase_plot(f, B, g₁, h₁, quiver_scaling, Δt, X0, solver = DifferentialEquations.Tsit5())
    p = ImplicitPlots.implicit_plot(B*g₁*h₁; xlims=(-2, 3), ylims=(-2.5, 2.5), resolution = 1000, label="")
    Plots.plot(p)
    ∇(vx, vy) = [fi(x[1] => vx, x[2] => vy) for fi in f]
    ∇pt(v, p, t) = ∇(v[1], v[2])
    function traj(v0)
        tspan = (0.0, Δt)
        prob = DifferentialEquations.ODEProblem(∇pt, v0, tspan)
        return DifferentialEquations.solve(prob, solver, reltol=1e-8, abstol=1e-8)
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