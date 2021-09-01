using DynamicPolynomials
@polyvar x[1:3]

poly = x[1]^2 - 2x[1]*x[2] + 3x[2]^2 - 2x[1]^2*x[2] + 2x[1]^2*x[2]^2 - 2x[2]*x[3] + 6x[3]^2 + 18x[2]^2*x[3] - 54x[2]*x[3]^2 + 142x[2]^2*x[3]^2

import CSDP
solver = CSDP.Optimizer
using SumOfSquares
function sos_min(sparsity)
    model = Model(solver)
    @variable(model, t)
    @objective(model, Max, t)
    con_ref = @constraint(model, poly - t in SOSCone(), sparsity = sparsity)
    optimize!(model)
    return value(t), moment_matrix(con_ref)
end

bound, ν = sos_min(NoSparsity())
bound

extractatoms(ν, 1e-6)

ν.basis

bound, ν = sos_min(MonomialSparsity())
bound

[sub.basis for sub in ν.sub_moment_matrices]

bound, ν = sos_min(MonomialSparsity(ChordalCompletion()))
bound

[sub.basis for sub in ν.sub_moment_matrices]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

