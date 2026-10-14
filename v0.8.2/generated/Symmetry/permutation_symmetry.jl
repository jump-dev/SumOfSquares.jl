import MutableArithmetics as MA
using MultivariatePolynomials
using MultivariateBases

using DynamicPolynomials
@polyvar x[1:4]

poly = sum(x) + sum(x.^2)

using SumOfSquares

using PermutationGroups
G = PermGroup([perm"(1,2,3,4)"])

import Clarabel
solver = Clarabel.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, Symmetry.VariablePermutation())
con_ref = @constraint(model, poly - t in SOSCone(), symmetry = pattern)
optimize!(model)
value(t)

gram_matrix(con_ref)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
