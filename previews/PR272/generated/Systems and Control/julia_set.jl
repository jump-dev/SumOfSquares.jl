function julia_map(point, c)
    a, b = point
    return [a^2 - b^2 + real(c), 2a * b + imag(c)]
end

escape_radius(c) = (1 + √(1 + 4 * abs(c))) / 2

using LinearAlgebra
function in_set(x, c, m=2000)
    r = escape_radius(c)
    for i in 1:m
        if norm(x) > r
            return false
        end
        x = julia_map(x, c)
    end
    return true
end

using SpecialFunctions
using DynamicPolynomials
β(α) = (α + 1) / 2
function circle_integral(mono::AbstractMonomial)
    if any(isodd, exponents(mono))
        return 0.0
    else
        return 2 * prod(gamma ∘ β, exponents(mono)) / gamma(sum(β, exponents(mono)))
    end
end
function disk_integral(mono::AbstractMonomial, r)
    d = degree(mono) + nvariables(mono)
    return circle_integral(mono) * r^d / d
end
function disk_integral(p::AbstractPolynomialLike, r)
    return sum(MultivariatePolynomials.coefficient(t) * disk_integral(monomial(t), r) for t in terms(p))
end

using SumOfSquares
function outer_approximation(solver, d::Int, c; α = 1/2)
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
    return value(v), value(w)
end

using ImplicitPlots
using Plots
function julia_plot(poly, c, n=200, m=1000; tol=1e-6, res = 1000)
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

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

c = -0.7 + 0.2im
v, w = outer_approximation(solver, 2, c)

julia_plot(v, c)

v, w = outer_approximation(solver, 4, c)

julia_plot(v, c)

c = -0.9 + 0.2im
v, w = outer_approximation(solver, 2, c)

julia_plot(v, c)

v, w = outer_approximation(solver, 4, c)

julia_plot(v, c)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

