using SumOfSquares

using DynamicPolynomials
@ncpolyvar x y
p = (x + 1.0im * y) * (x - im * y)

import CSDP
model = Model(CSDP.Optimizer)
add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge)
add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge)
add_bridge(model, SumOfSquares.COI.Bridges.Variable.HermitianToSymmetricPSDBridge)
add_bridge(model, SumOfSquares.COI.Bridges.Constraint.SplitZeroBridge)
MOI.Bridges.add_bridge(backend(model).optimizer, PolyJuMP.ZeroPolynomialBridge{Complex{Float64}})
MOI.Bridges.add_bridge(backend(model).optimizer, SumOfSquares.Bridges.Constraint.SOSPolynomialBridge{Complex{Float64}})
cone = NonnegPolyInnerCone{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}()
certificate = SumOfSquares.Certificate.MaxDegree(cone, MonomialBasis, 2)
c = SumOfSquares.add_constraint(model, p, cone, ideal_certificate = certificate)
optimize!(model)
sos_decomposition(c, 1e-6) #!src

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

