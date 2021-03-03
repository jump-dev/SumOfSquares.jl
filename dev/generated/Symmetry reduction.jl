import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using DynamicPolynomials
@polyvar x[1:4]

poly = sum(x) + sum(x.^2)

using SumOfSquares
include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "symmetry.jl"))

using PermutationGroups
function action(mono::AbstractMonomial, p::Perm)
    v = variables(mono)
    MP.substitute(MP.Eval(), mono, v => [v[i^p] for i in eachindex(v)])
end
G = PermGroup([perm"(1,2,3,4)"])

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G, action))
optimize!(model)
value(t)

gram_matrix(con_ref)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

