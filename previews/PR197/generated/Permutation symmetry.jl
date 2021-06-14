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

using PermutationGroups
G = PermGroup([perm"(1,2,3,4)"])

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
certificate = Certificate.SymmetricIdeal(Certificate.MaxDegree(SOSCone(), MonomialBasis, maxdegree(poly)), G, Certificate.VariablePermutation())
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = certificate)
optimize!(model)
value(t)

for g in gram_matrix(con_ref).sub_gram_matrices
    println(g.basis.polynomials)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

