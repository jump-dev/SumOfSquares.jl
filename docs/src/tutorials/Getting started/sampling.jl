# # Sampling basis

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Getting started/sampling.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Getting started/sampling.ipynb)
# **Contributed by**: Benoît Legat

using Test #src
using DynamicPolynomials
using SumOfSquares
import MultivariateBases as MB
import SDPLR
import Hypatia
import SCS
import BMSOS

# In this tutorial, we show how to use a different polynomial basis
# for enforcing the equality between the polynomial and its Sum-of-Squares decomposition.

@polyvar x
p = x^4 - 4x^3 - 2x^2 + 12x + 3

scs = SCS.Optimizer
sdplr = optimizer_with_attributes(SDPLR.Optimizer, "maxrank" => (m, n) -> 4)
hypatia = Hypatia.Optimizer
bmsos = BMSOS.Optimizer
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
    if !feasible
        @test value(γ) ≈ -6 rtol=1e-4
    end
end
test(scs)
test(sdplr)
test(hypatia)
test(bmsos, true)

function test_rand(solver, d, B)
    model = Model(solver)
    set_silent(model)
    p = MB.algebra_element(rand(2d+1), MB.SubBasis{B}(monomials(x, 0:2d)))
    @constraint(model, p in SOSCone(), zero_basis = BoxSampling([-1], [1]))
    optimize!(model)
    return solve_time(model)
end

test_rand(scs, 2, MultivariateBases.Trigonometric)
test_rand(bmsos, 4, MultivariateBases.Trigonometric)
