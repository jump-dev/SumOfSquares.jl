# # Outer approximation of Julia set

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Systems and Control/julia_set.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Systems and Control/julia_set.ipynb)
# **Adapted from**: Section 7.1.3 of [KHJ14]
#
# [KHJ14] Milan Korda, Didier Henrion, and Colin N. Jones.
# *Convex computation of the maximum controlled invariant set for polynomial control systems*.
# SIAM Journal on Control and Optimization 52.5 (2014): 2944-2969.

# The Julia map is defined as follows:

function julia_map(point, c)
    a, b = point
    return [a^2 - b^2 + real(c), 2a * b + imag(c)]
end


# The *escape radius" is the radius `r` such that `r^2 ≥ r + abs(c)`.
# Ouside of the circle of that radius, all points diverge so we know the Julia set belongs to that circle.

escape_radius(c) = (1 + √(1 + 4 * abs(c))) / 2

# The following function implements [KHJ14, (8)].

using SumOfSquares
using DynamicPolynomials
function outer_approximation(solver, d::Int; c = -0.7 + 0.2im, α = 1/2)
    @polyvar x[1:2]
    model = SOSModel(solver)
    r = escape_radius(c)
    S = @set sum(x.^2) <= r^2
    @variable(model, v, Poly(monomials(x, 0:2d)))
    @variable(model, w0, SOSPoly(monomials(x, 0:d)))
    @variable(model, w1, SOSPoly(monomials(x, 0:(d - 1))))
    @constraint(model, α * v(x => julia_map(x, c)) <= v, domain = S)
    w = w0 + w1 * (r^2 - sum(x.^2))
    @constraint(model, w >= v + 1, domain = S)
    @objective(model, Min, disk_integral(w, r))
    optimize!(model)
    println(solution_summary(model))
    if primal_status(model) == MOI.NO_SOLUTION
        return
    end
    return value(v), value(w), model
end

# The following function plots the Julia set with the outer approximation.

using ImplicitPlots
using Plots
function julia_plot(poly, n=400, m=1000; c = DEFAULT_C, tol=1e-6, res = 1000)
    r = escape_radius(c)
    p = implicit_plot(poly; xlims=(-r, r) .* 1.1, ylims=(-r, r), resolution = res, label="")
    θ = range(0, stop=2π, length=100)
    points = Vector{Float64}[]
    as = range(-r, r, length=n)
    bs = range(-r, r, length=n)
    for a in as, b in bs
        point = [a, b]
        if in_set(point, c, m)
            push!(points, point)
        end
    end
    xs = [point[1] for point in points]
    ys = [point[2] for point in points]
    scatter!(p, xs, ys, label="", markerstrokewidth=0, markersize=1, m=:pixel)
    return p
end

# We need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.8/installation/#Supported-solvers) for a list of the available choices.

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

# Let's start with the value of `c` corresponding to the left image of [KHJ14, Figure 3] and with degree 2.

c = -0.7 + 0.2im
v, w, model = outer_approximation(solver, 2; c)
julia_plot(v)

# Let's now look at degree 4.

v, w, model = outer_approximation(solver, 4; c)
julia_plot(v)

# Let's finish with degree 6.

v, w, model = outer_approximation(solver, 6; c)
julia_plot(v)

# Let's now use the value of `c` corresponding to the right image of [KHJ14, Figure 3] and with degree 2.

c = -0.9 + 0.2im
v, w, model = outer_approximation(solver, 2; c)
julia_plot(v)

# Let's now look at degree 4.

v, w, model = outer_approximation(solver, 4; c)
julia_plot(v)

# Let's finish with degree 6.

v, w, model = outer_approximation(solver, 6; c)
julia_plot(v)
