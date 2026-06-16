# # Sampling basis

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/sampling.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/sampling.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
using DynamicPolynomials
using SumOfSquares
import MultivariateBases as MB

# In this tutorial, we show how to use a different polynomial basis
# for enforcing the equality between the polynomial and its Sum-of-Squares decomposition.

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 3

# We want to find the minimum of the above polynomial (which is -6).

model = Model()
set_silent(model)
γ = -6
@variable(model, γ)
@objective(model, Max, γ)
@constraint(model, p - γ in SOSCone(), zero_basis = BoxSampling([-1], [1]))
set_optimizer()
optimize!(model)

function test(solver, feas::Bool)
    model = Model(solver)
    set_silent(model)
    if feas
        γ = -6
    else
        @variable(model, γ)
        @objective(model, Max, γ)
    end
    @constraint(model, p - γ in SOSCone(), zero_basis = BoxSampling([-1], [1]))
    optimize!(model)
    @test primal_status == MOI.FEASIBLE_POINT
    if !feas
        @test value(γ) ≈ -6 rtol=1e-4
    end
end


import SDPLR

import BMSOS

import LowRankOpt as LRO
import Percival
import Dualization

sdplr = optimizer_with_attributes(SDPLR.Optimizer, "maxrank" => (m, n) -> 4)
bmsos = BMSOS.Optimizer
bmlbfgs = Dualization.dual_optimizer(
    optimizer_with_attributes(
        LRO.Optimizer,
        "solver" => LRO.BurerMonteiro.Solver,
        "sub_solver" => Percival.PercivalSolver,
        "ranks" => [4],
        "square_scalars" => true,
    );
    assume_min_if_feasibility = true,
)
test(sdplr)
test(bmsos, true)
test(bmlbfgs, true)

function bench_rand(solver, d, B)
    model = Model(solver)
    set_silent(model)
    p = MB.algebra_element(rand(2d+1), MB.SubBasis{B}(monomials(x, 0:2d)))
    @constraint(model, p in SOSCone(), zero_basis = BoxSampling([-1], [1]))
    optimize!(model)
    return solve_time(model)
end

import SCS
scs = SCS.Optimizer
bench_rand(scs, 100, MultivariateBases.Trigonometric)

import Hypatia
hypatia = Hypatia.Optimizer
test_rand(hypatia, 100, MultivariateBases.Trigonometric)

test_rand(sdplr, 100, MultivariateBases.Trigonometric)
test_rand(bmsos, 100, MultivariateBases.Trigonometric)
test_rand(bmlbfgs, 100, MultivariateBases.Trigonometric)
