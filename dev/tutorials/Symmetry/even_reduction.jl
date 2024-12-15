# # Even reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry/even_reduction.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry/even_reduction.ipynb)

using Test #src
using DynamicPolynomials
@polyvar x

# We would like to find the minimum value of the following polynomial:

poly = x^4 - 2x^2

using SumOfSquares

# We define the custom action as follows:

struct OnSign <: Symmetry.OnMonomials end
using PermutationGroups
import SymbolicWedderburn
SymbolicWedderburn.coeff_type(::OnSign) = Float64
function SymbolicWedderburn.action(::OnSign, p::Permutation, mono::AbstractMonomial)
    if isone(p) || iseven(DynamicPolynomials.degree(mono))
        return 1 * mono
    else
        @assert p.perm == perm"(1,2)"
        return -1 * mono
    end
end
G = PermGroup([perm"(1,2)"])

# We can exploit the symmetry as follows:

import Clarabel
solver = Clarabel.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
pattern = Symmetry.Pattern(G, OnSign())
con_ref = @constraint(model, poly - t in SOSCone(), symmetry = pattern)
optimize!(model)
@test value(t) ≈ -1 #src
value(t)

# We indeed find `-1`, let's verify that symmetry was exploited:

@test length(gram_matrix(con_ref).blocks) == 2 #src
polys = gram_matrix(con_ref).blocks[1].basis.bases[].elements #src
@test polys[1] ≈ 1 #src
@test polys[2] ≈ x^2 #src
polys = gram_matrix(con_ref).blocks[2].basis.bases[].elements #src
@test polys[] ≈ x #src
gram_matrix(con_ref)
