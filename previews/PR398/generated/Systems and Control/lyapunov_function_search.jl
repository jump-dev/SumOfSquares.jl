using DynamicPolynomials
@polyvar x[1:3]

f = [-x[1]^3 - x[1] * x[3]^2,
     -x[2] - x[1]^2 * x[2],
     -x[3] - 3x[3] / (x[3]^2 + 1) + 3x[1]^2 * x[3]]

using SumOfSquares
using CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

monos = x.^2

@variable(model, V, Poly(monos))

@constraint(model, V >= sum(x.^2))

using LinearAlgebra # Needed for `dot`
dVdt = dot(differentiate(V, x), f)

P = dVdt.num

@constraint(model, P <= 0)

JuMP.optimize!(model)

JuMP.primal_status(model)

value(V)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
