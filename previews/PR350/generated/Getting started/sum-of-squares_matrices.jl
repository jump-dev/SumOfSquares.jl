using DynamicPolynomials
@polyvar x
P = [x^2 - 2x + 2 x
            x     x^2]

using SumOfSquares

import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = SOSModel(solver)
mat_cref = @constraint(model, P in PSDCone())
optimize!(model)
termination_status(model)

@polyvar y[1:2]
p = vec(y)' * P * vec(y)

X = monomials(p)
unipartite = Certificate.NewtonDegreeBounds(tuple())

multipartite = Certificate.NewtonDegreeBounds(([x], y))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
