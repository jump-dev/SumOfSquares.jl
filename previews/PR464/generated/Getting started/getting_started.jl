using DynamicPolynomials
@polyvar x y
p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4

using SumOfSquares
import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver)
con_ref = @constraint(model, p >= 0)
optimize!(model)
primal_status(model)

q = gram_matrix(con_ref)

sosdec = SOSDecomposition(q)

p = 4*x^4*y^6 + x^2 - x*y^2 + y^2

model = Model(solver)
con_ref = @constraint(model, p in SOSCone())
optimize!(model)
primal_status(model)

sos_decomposition(con_ref)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
