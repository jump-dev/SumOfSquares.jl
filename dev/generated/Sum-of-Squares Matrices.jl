using DynamicPolynomials
@polyvar x
P = [x^2 - 2x + 2 x
            x     x^2]

using SumOfSquares

import CSDP
factory = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = SOSModel(factory)
mat_cref = @constraint(model, P in PSDCone())
optimize!(model)
termination_status(model) #!src

@polyvar y[1:2]
p = vec(y)' * P * vec(y)

X = monomials(p)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

