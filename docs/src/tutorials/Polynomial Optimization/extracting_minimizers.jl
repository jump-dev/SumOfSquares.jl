# # Extracting minimizers

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Polynomial Optimization/extracting_minimizers.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Polynomial Optimization/extracting_minimizers.ipynb)
# **Adapted from**: Example 6.23 of [L09]
#
# [L09] Laurent, Monique.
# *Sums of squares, moment matrices and optimization over polynomials.*
# Emerging applications of algebraic geometry (2009): 157-270.
# World Scientific, **2009**.

# ## Introduction

# Consider the polynomial optimization problem [L09, Example 6.23] of
# minimizing the linear function $-x_1 - x_2$
# over the basic semialgebraic set defined by the inequalities
# $x_2 \le 2x_1^4 - 8x_1^3 + 8x_1^2 + 2$,
# $x_2 \le 4x_1^4 - 32x_1^3 + 88x_1^2 - 96x_1 + 36$ and the box constraints
# $0 \le x_1 \le 3$ and $0 \le x_2 \le 4$,

using Test #src
using DynamicPolynomials
DynamicPolynomials.@polyvar x[1:2]
p = -sum(x)
using SumOfSquares
K = @set x[1] >= 0 && x[1] <= 3 && x[2] >= 0 && x[2] <= 4 && x[2] <= 2x[1]^4 - 8x[1]^3 + 8x[1]^2 + 2 && x[2] <= 4x[1]^4 - 32x[1]^3 + 88x[1]^2 - 96x[1] + 36
f1 = 2x[1]^4 - 8x[1]^3 + 8x[1]^2 + 2 
f2 = 4x[1]^4 - 32x[1]^3 + 88x[1]^2 - 96x[1] + 36

import MultivariatePolynomials as MP

function nonlinear_vars!(vars, p)
    for mono in MP.monomials(p)
        if MP.degree(mono) > 1
            for var in MP.effective_variables(mono)
                push!(vars, var)
            end
        end
    end
end

function nonlinear_vars(K::BasicSemialgebraicSet{T,P}) where {T,P}
    vars = Set{MP.variable_union_type(P)}()
    for p in inequalities(K)
        nonlinear_vars!(vars, p)
    end
    for p in equalities(K)
        nonlinear_vars!(vars, p)
    end
    return sort(collect(vars), rev=true)
end

function polyhedral_proj!(vars, v, K)
    nl_vars = nonlinear_vars(K)
    lin_vars = setdiff(vars, nl_vars)
end

function proj!(v)
    v[1] = min(max(v[1], 0), 3)
    v[2] = min(v[2], f1(v[1]))
    v[2] = min(v[2], f2(v[1]))
    v[2] = min(max(v[2], 0), 4)
    return v
end

proj!([3, 4])

function laurent(ν)
    μ = measure(ν)
    v = [moment_value(μ, var) for var in x]
    @show v
    proj!(v)
    return v, p(x => v)
end

import MultivariatePolynomials as MP
function gaussian(ν, leading::Bool = true)
    monos = MP.monomials(MP.variables(ν.basis.monomials), 0:1)
    I = [MultivariateMoments._index(ν.basis, mono) for mono in monos]
    Q = ν.Q[I, I]
    F = eigen(Q)
    display(F)
    best_v = nothing
    best_obj = Inf
    J = leading ? size(F.vectors, 2) : axes(F.vectors, 2)
    for i in J
        v = F.vectors[:, i]
        @show v
        v /= v[1]
        v = v[2:end]
        @show v
        proj!(v)
        @show v
        obj = p(x => v)
        @show obj
        if obj < best_obj
            best_v = v
            best_obj = obj
        end
    end
    return best_v, best_obj
end

# We will now see how to find the optimal solution using Sum of Squares Programming.
# We first need to pick an SDP solver, see [here](https://jump.dev/JuMP.jl/v1.12/installation/#Supported-solvers) for a list of the available choices.

import Clarabel
solver = Clarabel.Optimizer

# A Sum-of-Squares certificate that $p \ge \alpha$ over the domain `S`, ensures that $\alpha$ is a lower bound to the polynomial optimization problem.
# The following function searches for the largest lower bound and finds zero using the `d`th level of the hierarchy`.

function solve(d)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = K, maxdegree = d)
    optimize!(model)
    println(solution_summary(model))
    return model
end

# The first level of the hierarchy gives a lower bound of `-7``

model4 = solve(4)
@test objective_value(model4) ≈ -7 rtol=1e-4 #src
@test termination_status(model4) == MOI.OPTIMAL #src
ν4 = moment_matrix(model4[:c])
laurent(ν4)

import MultivariateMoments as MM
function _vec(Q)
    n = size(Q, 2)
    return [Q[i, j] for j in 1:n for i in 1:j]
end
function truncate(ν, d)
    vars = MP.variables(ν.basis)
    monos = MP.monomials(vars, 0:d)
    I = [MM._index(ν.basis, mono) for mono in monos]
    return MM.MomentMatrix(
        MM.SymMatrix(_vec(ν.Q[I, I]), length(I)),
        MonomialBasis(monos),
    )
end

# The second level improves the lower bound

model5 = solve(5)
@test objective_value(model5) ≈ -20/3 rtol=1e-4 #src
@test termination_status(model5) == MOI.OPTIMAL #src
ν5 = moment_matrix(model5[:c])
laurent(ν5)

# The third level finds the optimal objective value as lower bound...

model7 = solve(7)
@test objective_value(model7) ≈ -5.5080 rtol=1e-4 #src
@test termination_status(model7) == MOI.OPTIMAL #src

# ...and proves it by exhibiting the minimizer.

ν7 = moment_matrix(model7[:c])
laurent(ν7)
η = atomic_measure(ν7, 1e-3)
@test length(η.atoms) == 1 #src
@test η.atoms[1].center ≈ [2.3295, 3.1785] rtol=1e-4 #src

# We can indeed verify that the objective value at `x_opt` is equal to the lower bound.

opt = [2.3295, 3.1785]
x_opt = η.atoms[1].center
@test x_opt ≈ [2.3295, 3.1785] rtol=1e-4 #src
p(x_opt)
#x_opt = 

#atomic_measure(ν4, UserRank())

compute_support!(ν5, FixedRank(1))
ν5.support

ν5_2 = truncate(ν5, 1)
compute_support!(ν5_2, FixedRank(1))
ν5_2.support

[p(x => x_opt) for p in ν5.support.I.p]
[p(x => [2.6666, 1.23456]) for p in ν5.support.I.p]
compute_support!(ν5, FixedRank(2))
ν5.support
pp = -6.945530612461534 - 2.2636172932906113*x[2] + x[2]^2

lag, system = PolyJuMP.lagrangian_kkt(MOI.MIN_SENSE, 1.0 * p, K)

using Macaulay
νmax3 = moment_matrix(system.I.p, Clarabel.Optimizer, 3)
laurent(νmax3)
νmax4 = moment_matrix(system.I.p, Clarabel.Optimizer, 4)
laurent(νmax4)
