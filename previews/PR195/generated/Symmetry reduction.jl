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
function action(mono::AbstractMonomial, p::Permutation)
    v = variables(mono)
    MP.substitute(MP.Eval(), mono, v => [v[i^p] for i in eachindex(v)])
end
function action(term::MP.AbstractTerm, el::Permutation)
    return MP.coefficient(term) * action(MP.monomial(term), el)
end
function action(poly::MP.AbstractPolynomial, el::Permutation)
    return MP.polynomial([action(term, el) for term in MP.terms(poly)])
end

G = PermGroup([perm"(1,2,3,4)"])

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
certificate = SymmetricIdeal(Certificate.MaxDegree(SOSCone(), MonomialBasis, maxdegree(poly)), G, action)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = certificate)
optimize!(model)
value(t)

for g in gram_matrix(con_ref).sub_gram_matrices
    println(g.basis.polynomials)
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

