A = [
     0  0  1
     0 -1  0
    -2  1 -1
]
bz = [3, 0, -4] - [0, -1, -6]
y = [1.5, -0.5, -5]

using DynamicPolynomials
@polyvar x[1:3]
p = -2x[1] + x[2] - x[3]
using SumOfSquares
e = A * x - y
f = e'e - bz'bz / 4
K = @set sum(x) <= 4 && 3x[2] + x[3] <= 6 && f >= 0 && 0 <= x[1] && x[1] <= 2 && 0 <= x[2] && 0 <= x[3] && x[3] <= 3

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

model4 = solve(4)
nothing # hide

model6 = solve(6)
nothing # hide

model8 = solve(8)
nothing # hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
