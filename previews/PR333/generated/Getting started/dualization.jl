using DynamicPolynomials
using SumOfSquares

import SCS
@polyvar x
p = (x + 1)^2 * (x + 2)^2
model_scs = Model(SCS.Optimizer)
con_ref = @constraint(model_scs, p in SOSCone())
optimize!(model_scs)

using Dualization
model_dual_scs = Model(dual_optimizer(SCS.Optimizer))
@objective(model_dual_scs, Max, 0.0)
con_ref = @constraint(model_dual_scs, p in SOSCone())
optimize!(model_dual_scs)

print_active_bridges(model_scs)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
