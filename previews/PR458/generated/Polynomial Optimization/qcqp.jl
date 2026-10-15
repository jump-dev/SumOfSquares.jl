using DynamicPolynomials
@polyvar x[1:6]
centers = [2, 2, 1, 4, 1, 4]
weights = [25, 1, 1, 1, 1, 1]
p = -weights' * (x .- centers).^2
using SumOfSquares
K = @set x[1] >= 0 && x[2] >= 0 &&
    x[3] >= 1 && x[3] <= 5 &&
    x[4] >= 0 && x[4] <= 6 &&
    x[5] >= 1 && x[5] <= 5 &&
    x[6] >= 0 && x[6] <= 10 &&
    (x[3] - 3)^2 + x[4] >= 4 &&
    (x[5] - 3)^2 + x[6] >= 4 &&
    x[1] - 3x[2] <= 2 &&
    -x[1] + x[2] <= 2 &&
    x[1] + x[2] <= 6 &&
    x[1] + x[2] >= 2

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

model3 = solve(4)
nothing # hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
