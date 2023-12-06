using DynamicPolynomials
@polyvar x[1:8]
p = sum(x[1:3])
using SumOfSquares
K = @set 0.0025 * (x[4] + x[6]) <= 1 &&
    0.0025 * (-x[4] + x[5] + x[7]) <= 1 &&
    0.01 * (-x[5] + x[8]) <= 1 &&
    100x[1] - x[1] * x[6] + 8333.33252x[4] <= 250000/3 &&
    x[2] * x[4] - x[2] * x[7] - 1250x[4] + 1250x[5] <= 0 &&
    x[3] * x[5] - x[3] * x[8] - 2500x[5] + 1250000 <= 0 &&
    100 <= x[1] && x[1] <= 10000 &&
    1000 <= x[2] && x[2] <= 10000 &&
    1000 <= x[3] && x[3] <= 10000 &&
    10 <= x[4] && x[4] <= 1000 &&
    10 <= x[5] && x[5] <= 1000 &&
    10 <= x[6] && x[6] <= 1000 &&
    10 <= x[7] && x[7] <= 1000 &&
    10 <= x[8] && x[8] <= 1000

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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
