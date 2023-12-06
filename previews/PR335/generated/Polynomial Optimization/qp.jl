using LinearAlgebra
c = [42, 44, 45, 47, 47.5]
Q = 100I

using DynamicPolynomials
@polyvar x[1:5]
p = c'x - x' * Q * x / 2
using SumOfSquares
K = @set x[1] >= 0 && x[1] <= 1 &&
         x[2] >= 0 && x[2] <= 1 &&
         x[3] >= 0 && x[3] <= 1 &&
         x[4] >= 0 && x[4] <= 1 &&
         x[5] >= 0 && x[5] <= 1 &&
         10x[1] + 12x[2] + 11x[3] + 7x[4] + 4x[5] <= 40

import Clarabel
solver = Clarabel.Optimizer

function solve(d)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = K, maxdegree = d)
    optimize!(model)
    println(solution_summary(model))
    return model
end

model2 = solve(2)
nothing # hide

model3 = solve(3)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
