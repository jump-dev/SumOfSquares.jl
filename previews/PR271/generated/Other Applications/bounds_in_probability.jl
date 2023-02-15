μ0 = 1

μ1  = 1

σ = 1/2

μ2 = σ^2 + μ1^2

using DynamicPolynomials
@polyvar x
monos = [1, x, x^2]
using SumOfSquares
μ = measure([μ0, μ1, μ2], monos)

using CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver);

@variable(model, poly, Poly(monos))

K = @set 0 <= x && x <= 5
con_ref = @constraint(model, poly >= 0, domain = K)

@constraint(model, poly >= 1, domain = (@set 4 <= x && x <= 5))

using LinearAlgebra
@objective(model, Min, poly ⋅ μ)

optimize!(model)
primal_status(model)

objective_value(model)

value(poly) * 37^2

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

