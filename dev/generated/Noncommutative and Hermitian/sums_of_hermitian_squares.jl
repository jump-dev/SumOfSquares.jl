using SumOfSquares

using DynamicPolynomials
@ncpolyvar x y
p = (x + 1.0im * y) * (x - im * y)

import CSDP
model = Model(CSDP.Optimizer)
MOI.Bridges.add_bridge(backend(model).optimizer, PolyJuMP.Bridges.Constraint.ZeroPolynomialBridge{Complex{Float64}})
MOI.Bridges.add_bridge(backend(model).optimizer, SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Complex{Float64}})
cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
con_ref = @constraint(model, p in cone)
optimize!(model)
sos_decomposition(con_ref, 1e-6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

