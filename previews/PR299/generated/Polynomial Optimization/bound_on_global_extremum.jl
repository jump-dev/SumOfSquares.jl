using DynamicPolynomials
@polyvar x1 x2

f1 = x1 + x2 + 1
f2 = 19 - 14x1 + 3x1^2 - 14x2 + 6x1*x2 + 3x2^2
f3 = 2x1 - 3x2
f4 = 18 - 32x1 + 12x1^2 + 48x2 - 36x1*x2 + 27x2^2
f = (1 + f1^2 * f2) * (30 + f3^2 * f4)

using SumOfSquares
using CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

@variable(model, γ)
@objective(model, Max, γ)

@constraint(model, f >= γ)

JuMP.optimize!(model)

JuMP.primal_status(model)

objective_value(model)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

