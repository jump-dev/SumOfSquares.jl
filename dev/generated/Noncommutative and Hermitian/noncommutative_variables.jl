using DynamicPolynomials
@ncpolyvar x y
p = (x * y + x^2)^2

using SumOfSquares
import CSDP
optimizer_constructor = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = Model(optimizer_constructor)
con_ref = @constraint(model, p in SOSCone())

optimize!(model)

certificate_basis(con_ref)

gram_matrix(con_ref).Q

sos_decomposition(con_ref) #!src

sos_decomposition(con_ref, 1e-6)       #!src

using DynamicPolynomials
@ncpolyvar x y
n = 10
p = (x + x^n * y^(2n) * x^n)^2

using SumOfSquares
model = Model(optimizer_constructor)
con_ref = @constraint(model, p in SOSCone())

optimize!(model)

certificate_basis(con_ref)

gram_matrix(con_ref).Q

sos_decomposition(con_ref, 1e-6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

