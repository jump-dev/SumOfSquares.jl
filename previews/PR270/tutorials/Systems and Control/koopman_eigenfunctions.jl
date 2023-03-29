# # Eigenfunctions of the Koopman operator

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Systems and Control/koopman_eigenfunctions.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Systems and Control/koopman_eigenfunctions.ipynb)
# **Adapted from**: Example 2.6 of [MAI20]
#
# [MAI20] Mauroy, Alexandre, Aivar Sootla, and Igor Mezić.
# *Koopman framework for global stability analysis.*
# The Koopman Operator in Systems and Control: Concepts, Methodologies, and Applications (2020): 35-58.

using Test #src
using DynamicPolynomials
@polyvar x[1:2]

a = 1
I = 0.05
ε = 0.08
γ = 1
F0 = [-x[2] - x[1] * (x[1] - 1) * (x[1] - a) + I, ε * (x[1] - γ * x[2])]

# We move equilibrium `(0.0256, 0.0256)` to the origin

x1 = x2 = 0.0256
F = [f(x => [x[1] + x1, x[2] + x2]) for f in F0]

# We compute the Jacobian at the equilibrium

J = [j(x => zeros(2)) for j in differentiate(F, x)]

# We see that its eigenvalues indeed have negative real part:
using LinearAlgebra
E = eigen(J)

# We set `w` as its dominant eigenvector:

λ = E.values[end]
w = E.vectors[:, end]

using SumOfSquares
r = 0.3
X = @set x[1]^2 + x[2]^2 ≤ r^2

# We define the the program for the FitzHugh-Nagumo problem below:
# `N` is the degree of `ϕN` as defined in [MAI20, p. 50] and `M` is the
# degree of the multipliers for the constraints of `X`.
# As `maxdegree` corresponds to the degree of the multiplier multiplied
# by `r^2 - x[1]^2 - x[2]^2`, we set it to `M + 2`.

function fitzhugh_nagumo(solver, N, M)
    model = SOSModel(solver)
    @variable(model, γ)
    @objective(model, Min, γ)
    @variable(model, ϕN, Poly(monomials(x, 2:N)))
    ϕ = w ⋅ x + ϕN
    ∇ϕ = differentiate(ϕ, x)
    @constraint(model, -γ ≤ F ⋅ ∇ϕ - λ * ϕ, domain = X, maxdegree = M + 2)
    @constraint(model, F ⋅ ∇ϕ - λ * ϕ ≤ γ, domain = X, maxdegree = M + 2)
    optimize!(model)
    return ϕ, model
end

# In [MAI20, p. 50], we read that the result is obtained with `N = 10` and
# `M = 20`.

import CSDP
ϕ, model = fitzhugh_nagumo(CSDP.Optimizer, 10, 20)
solution_summary(model)

# The optimal value of `ϕ` is obtained as follows:

ϕ_opt = value(ϕ)

# Its plot is given below:

using Plots
x1s = x2s = range(-0.3, stop = 0.3, length = 40)
ϕs = ϕ_opt.(x1s', x2s)
contour(x1s, x2s, ϕs, levels=[-0.2, -0.1, 0, 0.1, 0.2, 0.5, 0.7])
θ = range(0, stop = 2π, length = 100)
plot!(r * cos.(θ), r * sin.(θ), label = "")
scatter!([0], [0], label = "")

# We can compute the Laplace average as follows:

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

# The error is given by:

norm(laplace - ϕs)