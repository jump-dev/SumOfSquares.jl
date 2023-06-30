using SumOfSquares

using DynamicPolynomials
@ncpolyvar x y
p = (x + im * y) * (x - im * y)

import CSDP
model = Model(CSDP.Optimizer)
cone = NonnegPolyInnerCone{MOI.HermitianPositiveSemidefiniteConeTriangle}()
con_ref = @constraint(model, p in cone)
optimize!(model)
sos_decomposition(con_ref, 1e-6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

