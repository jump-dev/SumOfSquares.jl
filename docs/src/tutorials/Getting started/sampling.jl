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
import DSDP
dsdp = DSDP.Optimizer
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
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @show objective_value(model)
    @show value(γ)
    if !feas
        if !isapprox(value(γ), -6, rtol=1e-3)
            @warn("$(value(γ)) != -6")
        end
        #@test value(γ) ≈ -6 rtol=1e-4
    end
    return model
end
model = test(scs, false);
model = test(sdplr, false);
test(hypatia, false);
#test(bmsos, true)
model = test(dsdp, false);

import Random
import TrigPolys
# See https://codeocean.com/capsule/8311748/tree/v1
function random_positive_poly(n; tol=1e-5)
    Random.seed!(0)
    p = TrigPolys.random_trig_poly(n)
    #N = 10000000
    N = 1000000
    p - minimum(TrigPolys.evaluate(TrigPolys.pad_to(p, N))) + n * tol
    a = zeros(2n + 1)
    a[1] = p.a0
    a[2:2:2n] = p.ac
    a[3:2:(2n+1)] = p.as
    return MB.algebra_element(
        a,
        MB.SubBasis{MB.Trigonometric}(monomials(x, 0:2n)),
    )
end
random_positive_poly(20)

function test_rand(solver, d, B)
    model = Model(solver)
    set_silent(model)
    p = if B == MB.Trigonometric
        random_positive_poly(d)
    else
        MB.algebra_element(rand(2d+1), MB.SubBasis{B}(monomials(x, 0:2d)))
    end
    @constraint(model, p in SOSCone(), zero_basis = BoxSampling([-1], [1]))
    optimize!(model)
    return solve_time(model)
end

using DataFrames
df = DataFrame(solver=String[], degree=Int[], time=Float64[])

function timing(solver, d, dual::Bool = false)
    name = MOI.get(MOI.instantiate(solver), MOI.SolverName())
    if dual
        solver = dual_optimizer(solver)
    end
    time = test_rand(solver, d, MB.Trigonometric)
    push!(df, (name, d, time))
end

using Dualization
using MosekTools
mosek = Mosek.Optimizer

d = 400
timing(bmsos, d)
timing(hypatia, d)
timing(dsdp, d)
timing(mosek, d)
#timing(sdplr, d)
#timing(dual_optimizer(scs), d)
#timing(scs, d, true)
#timing(mosek, d, true)

using Printf

function table(degs, solvers)
    print("|     |")
    for solver in solvers
        print(" $solver |")
    end
    println()
    for _ in 0:length(solvers)
        print("|-----")
    end
    println("|")
    for deg in degs
        print("| ", deg, " |")
        for solver in solvers
            times = subset(df, :solver => c -> c .== Ref(solver), :degree => d -> d .== deg).time
            if isempty(times)
                print("  |")
            else
                @printf(" %.3e |", minimum(times))
            end
        end
        println()
    end
end

table([100, 200, 300, 400, 500, 800], ["BMSOS", "Hypatia", "DSDP", "Mosek"])
