using DynamicPolynomials
@polyvar x[1:3]

poly = 1 + x[1]^4 + x[1] * x[2] + x[2]^4 + x[3]^2

import CSDP
solver = CSDP.Optimizer
using SumOfSquares
function sos_check(sparsity)
    model = Model(solver)
    con_ref = @constraint(model, poly in SOSCone(), sparsity = sparsity)
    optimize!(model)
    println(solution_summary(model))
    return gram_matrix(con_ref)
end

g = sos_check(Sparsity.NoPattern())
g.basis.monomials

g = sos_check(Sparsity.SignSymmetry())
monos = [sub.basis.monomials for sub in g.blocks]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

