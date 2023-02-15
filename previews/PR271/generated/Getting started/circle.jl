using DynamicPolynomials
using SumOfSquares
@polyvar x y
S = @set x^2 + y^2 == 1

import CSDP
model = SOSModel(CSDP.Optimizer)
set_silent(model)
con_ref = @constraint(model, 1 - y^2 >= 0, domain = S)
optimize!(model)

solution_summary(model)

sos_decomposition(con_ref, 1e-6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

