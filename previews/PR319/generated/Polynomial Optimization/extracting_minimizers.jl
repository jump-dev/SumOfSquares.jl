using DynamicPolynomials
@polyvar x[1:2]
p = -sum(x)
using SumOfSquares
K = @set x[1] >= 0 && x[1] <= 3 && x[2] >= 0 && x[2] <= 4 && x[2] <= 2x[1]^4 - 8x[1]^3 + 8x[1]^2 + 2 && x[2] <= 4x[1]^4 - 32x[1]^3 + 88x[1]^2 - 96x[1] + 36

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

function solve(d)
    model = SOSModel(solver)
    @variable(model, α)
    @objective(model, Max, α)
    @constraint(model, c, p >= α, domain = K, maxdegree = d)
    optimize!(model)
    println(solution_summary(model))
    return model
end

model4 = solve(4)

model5 = solve(5)

model7 = solve(7)

ν7 = moment_matrix(model7[:c])
η = atomic_measure(ν7, 1e-3) # Returns nothing as the dual is not atomic

x_opt = η.atoms[1].center
p(x_opt)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
