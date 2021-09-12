using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
using SumOfSquares
S = @set x >= 0 && y >= 0 && x + y >= 1
p(x=>1, y=>0), p(x=>1//2, y=>1//2), p(x=>0, y=>1)

import Ipopt
model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
optimize!(model)
@show termination_status(model)
@show value(a)
@show value(b)
@show objective_value(model)

using Ipopt
model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
peval(a, b) = p(x=>a, y=>b)
register(model, :peval, 2, peval, autodiff=true)
@NLobjective(model, Min, peval(a, b))
optimize!(model)
@show termination_status(model)
@show value(a)
@show value(b)
@show objective_value(model)

# Sum-of-Squares approach

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c3, p >= α, domain = S)
optimize!(model)
@show termination_status(model)
@show objective_value(model)

ν3 = moment_matrix(c3)
extractatoms(ν3, 1e-3) # Returns nothing as the dual is not atomic

model = SOSModel(solver)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c5, p >= α, domain = S, maxdegree = 5)
optimize!(model)
@show termination_status(model)
@show objective_value(model)

ν5 = moment_matrix(c5)
extractatoms(ν5, 1e-3)

ν3 = moment_matrix(c3)
SumOfSquares.MultivariateMoments.computesupport!(ν3, 1e-3)

ν5 = moment_matrix(c5)
SumOfSquares.MultivariateMoments.computesupport!(ν5, 1e-3)

SemialgebraicSets.computegröbnerbasis!(ideal(ν5.support))
ν5.support

using HomotopyContinuation
solver = SemialgebraicSetsHCSolver(; excess_residual_tol = 2e-2, real_tol = 2e-2, compile = false)
extractatoms(ν5, 1e-3, solver)

F = HomotopyContinuation.System(ν5.support)
res = HomotopyContinuation.solve(F, solver.options...)
path_results(res)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

