# # Eigenfunctions of the Koopman operator

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Systems and Control/koopman_eigenfunctions.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Systems and Control/koopman_eigenfunctions.ipynb)
# **Adapted from**: Example 2.6 of [MAI20]
#
# [MAI20] Mauroy, Alexandre, Aivar Sootla, and Igor Mezić.
# *Koopman framework for global stability analysis.*
# The Koopman Operator in Systems and Control: Concepts, Methodologies,
# and Applications (2020): 35-58.

using Test #src
using DynamicPolynomials
using LinearAlgebra
@polyvar x[1:2]

a = 1
I = 0.05
ε = 0.08
γ = 1
F0 = [-x[2] - x[1] * (x[1] - 1) * (x[1] - a) + I, ε * (x[1] - γ * x[2])]

# We move the equilibrium, approximately `(0.0256, 0.0256)`, to the origin.
# Since `γ = 1`, both coordinates are equal. Newton's method refines the
# rounded value so that the translated vector field vanishes at the origin.

x_eq = let z = 0.0256
    for _ in 1:5
        z -= (-z - z * (z - 1)^2 + I) / (-2 + 4z - 3z^2)
    end
    z
end
F = [f(x => x .+ x_eq) for f in F0]
@test norm([f(x => zeros(2)) for f in F]) < 1e-14 #src

# We compute the Jacobian at the equilibrium.

J = [j(x => zeros(2)) for j in differentiate(F, x)]

# Its eigenvalues have negative real part. The linear term `w ⋅ x` of a
# Koopman eigenfunction must satisfy `J' * w = λ * w`, so we need a **left**
# eigenvector of `J`. A right eigenvector leaves a nonzero linear residual
# which the higher-degree terms cannot cancel. Example 2.6 calls `w` a
# right eigenvector, but this contradicts Eq. (2.7) on p. 42 of the book.

E = eigen(J')
λ = E.values[end]
w = E.vectors[:, end]
@test all(<(0), E.values) #src
@test J' * w ≈ λ * w #src

# Eigenfunctions are defined up to a scalar factor. We choose a unit gradient
# with positive second component to fix the sign of the contour labels.
# Figure 2.3 appears to use the original coordinates, with the equilibrium
# at `(x_eq, x_eq)`, and a different normalization. Our axes are shifted to
# the equilibrium and our gradient has unit norm, so the contour labels
# need not coincide with those in the book.

if w[2] < 0
    w = -w
end

# We approximate the eigenfunction on the disk of radius `r`. Working in
# coordinates `y = x / r` avoids the very small high-degree monomials on
# this disk. The scaled vector field is `G(y) = F(r * y) / r`.

using SumOfSquares
r = 0.3
@polyvar y[1:2]
G = [f(x => r .* y) / r for f in F]
X = @set y[1]^2 + y[2]^2 ≤ 1

# We define the program for the FitzHugh-Nagumo problem below:
# `N` is the degree of `ϕN` as defined in [MAI20, p. 50] and `M` is the
# degree of the multipliers for the constraints of `X`.
# As `maxdegree` corresponds to the degree of the multiplier multiplied
# by `1 - y[1]^2 - y[2]^2`, we set it to `M + 2`.

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

# In [MAI20, p. 50], the result is obtained with `N = 10` and `M = 20`.

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
ϕ, model = fitzhugh_nagumo(solver, 10, 20)
solution_summary(model)
@test termination_status(model) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL] #src
@test primal_status(model) in [MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT] #src
@test 0 ≤ objective_value(model) < 1e-5 #src

# The optimal polynomial is obtained as follows:

ϕ_opt = value(ϕ)
residual_opt = F ⋅ differentiate(ϕ_opt, x) - λ * ϕ_opt
@test abs(ϕ_opt(x => zeros(2))) < 1e-14 #src
@test [d(x => zeros(2)) for d in differentiate(ϕ_opt, x)] ≈ w #src

# We plot its level sets inside the disk where the residual is constrained.
# Values outside this disk are extrapolations and are excluded from both
# the plot and the error calculation below.

using Plots
x1s = x2s = range(-r, r, length = 40)
inside = [u^2 + v^2 ≤ r^2 for v in x2s, u in x1s]
ϕs = [
    inside[i, j] ? ϕ_opt(x => [u, v]) : NaN for
    (i, v) in enumerate(x2s), (j, u) in enumerate(x1s)
]
residuals = [residual_opt(x => [u, v]) for v in x2s, u in x1s] #src
@test maximum(abs, residuals[inside]) ≤ objective_value(model) + 1e-7 #src
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

# We can check the Koopman identity using the finite-time Laplace average
# of `ϕ_opt`. An exact eigenfunction satisfies
# `exp(-λ * t) * ϕ(S(t, x)) = ϕ(x)` along a trajectory.
# Solve each trajectory once, then integrate using its dense interpolant.
# Tight tolerances matter because `exp(-λ * t)` amplifies integration errors.

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

# We report the maximum absolute discrepancy on the disk. Unlike a matrix
# norm over the enclosing square, this error excludes unconstrained points
# and does not grow merely because the grid is refined. It checks the
# eigenfunction identity over a finite time, rather than providing a bound
# on the error with respect to the exact eigenfunction.

maximum(abs, laplace[inside] - ϕs[inside])

# The disk is not forward invariant: some trajectories starting inside it
# leave it, so a small residual on the disk does not guarantee a small
# Laplace-average discrepancy for all those initial states. Increasing `T`
# does not fix extrapolation outside the approximation domain. For initial
# states in the smaller disk of radius `r / 2`, the discrepancy is much
# smaller:

inner = [u^2 + v^2 ≤ (r / 2)^2 for v in x2s, u in x1s]
maximum(abs, laplace[inner] - ϕs[inner])
@test maximum(abs, laplace[inner] - ϕs[inner]) < 1e-4 #src

# The book reports an overall error of `2.5e-3`, but does not specify the
# observable used for the Laplace average, its finite integration horizon,
# or the error norm. That number cannot be directly compared with the
# finite-time identity discrepancies above, especially when eigenfunctions
# have different scalar normalizations.
