using DynamicPolynomials
@polyvar x y
p = x^3 - x^2 + 2x*y -y^2 + y^3
using SumOfSquares
S = @set x >= 0 && y >= 0 && x + y >= 1
p(x=>1, y=>0), p(x=>1//2, y=>1//2), p(x=>0, y=>1)

import Ipopt
model = Model(Ipopt.Optimizer)
@variable(model, a >= 0)
@variable(model, b >= 0)
@constraint(model, a + b >= 1)
@NLobjective(model, Min, a^3 - a^2 + 2a*b - b^2 + b^3)
optimize!(model)

solution_summary(model)

value(a), value(b)

f(a, b) = p(x => a, y => b)
∇p = differentiate(p, [x, y])
function ∇f(g, a, b)
    for i in eachindex(g)
        g[i] = ∇p[i](x => a, y => b)
    end
end
∇²p = differentiate(∇p, [x, y])
function ∇²f(H, a, b)
    for j in axes(∇²p, 2)
        for i in j:size(∇²p, 1)
            H[i, j] = ∇²p[i, j](x => a, y => b)
        end
    end
end
using Ipopt
gmodel = Model(Ipopt.Optimizer)
@variable(gmodel, a >= 0)
@variable(gmodel, b >= 0)
@constraint(gmodel, a + b >= 1)
register(gmodel, :f, 2, f, ∇f, ∇²f)
@NLobjective(gmodel, Min, f(a, b))
optimize!(gmodel)

solution_summary(gmodel)

value(a), value(b)

import SCS
scs = SCS.Optimizer
import Dualization
dual_scs = Dualization.dual_optimizer(scs)

model = SOSModel(dual_scs)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c3, p >= α, domain = S)
optimize!(model)

solution_summary(model)

ν3 = moment_matrix(c3)
atomic_measure(ν3, 1e-3) # Returns nothing as the dual is not atomic

model = SOSModel(dual_scs)
@variable(model, α)
@objective(model, Max, α)
@constraint(model, c4, p >= α, domain = S, maxdegree = 4)
optimize!(model)

moment_matrix(c4)

function sos(solver, deg)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = S, maxdegree = deg, newton_polytope = nothing)
    optimize!(model)
    return model
end
dual_model4 = sos(dual_scs, 4)
nothing #hide

solution_summary(dual_model4)

dual_ν4 = moment_matrix(dual_model4[:c])

using LinearAlgebra
svdvals(Matrix(dual_ν4.Q))

atomic_measure(dual_ν4, FixedRank(4))

model4 = sos(scs, 4)
nothing #hide

solution_summary(model4)

ν4 = moment_matrix(model4[:c])

svdvals(Matrix(ν4.Q))

atomic_measure(ν4, FixedRank(3))

ν3 = moment_matrix(c3)
SumOfSquares.MultivariateMoments.compute_support!(ν3, LeadingRelativeRankTol(1e-3))

ν4 = moment_matrix(model4[:c])
SumOfSquares.MultivariateMoments.compute_support!(ν4, FixedRank(3))

SemialgebraicSets.compute_gröbner_basis!(ideal(ν4.support))
ν4.support

collect(ν4.support)

using HomotopyContinuation
algebraic_solver = SemialgebraicSetsHCSolver(; excess_residual_tol = 1e-1, real_tol = 1e-1, compile = false)
atomic_measure(ν4, FixedRank(3), Echelon(), algebraic_solver)

F = HomotopyContinuation.System(ν4.support)
res = HomotopyContinuation.solve(F, algebraic_solver.options...)
path_results(res)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
