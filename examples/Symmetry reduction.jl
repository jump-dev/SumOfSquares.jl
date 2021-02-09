# # Symmetry reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry reduction.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry reduction.ipynb)
# **Adapted from**: https://github.com/kalmarek/SymbolicWedderburn.jl/blob/tw/ex_sos/examples/ex_C4.jl

using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl"

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using Test #src
using DynamicPolynomials
@polyvar x[1:4]

# We would like to find the minimum value of the polynomial

poly = sum(x) + sum(x.^2)

# As we can decouple the problem for each `x[i]` for which `x[i] + x[i]^2` has
# minimum value 0.25, we would expect to get `-1` as answer.
# Can this decoupling be exploited by SumOfSquares as well ?
# For this, we need to use a certificate that can exploit the permutation symmetry of the polynomial.
# This is still a work in progress in SumOfSquares, so we need to define things here:

using SymbolicWedderburn
using PermutationGroups
using Cyclotomics
using SumOfSquares

function SymbolicWedderburn.ExtensionHomomorphism(basis::MB.MonomialBasis)
    monos = basis.monomials
    basis_exps = Vector{Vector{Int}}(undef, length(monos))
    basis_dict = Dict{Vector{Int}, Int}()
    sizehint!(basis_dict, length(monos))

    for (i, b) in enumerate(monos)
        e = MP.exponents(b) # so that we allocate exponents only once
        basis_exps[i] = e
        basis_dict[e] = i
    end

    return SymbolicWedderburn.ExtensionHomomorphism(basis_exps, basis_dict)
end
function MP.polynomialtype(::Type{<:MB.AbstractPolynomialVectorBasis{PT}}, T::Type) where PT
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end
struct SymmetricIdeal{CT, GT} <: Certificate.AbstractIdealCertificate
    cone::CT
    group::GT
end
SumOfSquares.matrix_cone_type(::Type{<:SymmetricIdeal{CT}}) where {CT} = SumOfSquares.matrix_cone_type(CT)
Certificate.get(::Type{<:SymmetricIdeal}, ::SumOfSquares.Certificate.GramBasisType) = Vector{MB.FixedPolynomialBasis}
Certificate.zero_basis_type(::Type{<:SymmetricIdeal}) = MB.MonomialBasis
Certificate.zero_basis(::SymmetricIdeal) = MB.MonomialBasis
Certificate.get(::SymmetricIdeal, ::Certificate.ReducedPolynomial, poly, domain) = poly
function Certificate.get(cert::SymmetricIdeal, ::Certificate.GramBasis, poly)
    basis = Certificate.maxdegree_gram_basis(MB.MonomialBasis, MP.variables(poly), MP.maxdegree(poly))
    R = SymbolicWedderburn.symmetry_adapted_basis(cert.group, basis)
    return map(R) do Ri
        FixedPolynomialBasis(convert(Matrix{Float64}, Ri) * basis.monomials)
    end
end
G = PermGroup([perm"(1,2,3,4)"])

# We can use this certificate as follows:

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G))
optimize!(model)
@test value(t) ≈ -1 #src
value(t)

# We indeed find `-1`, let's verify that symmetry was exploited:

@test length(gram_matrix(con_ref).sub_gram_matrices) == 3 #src
gram_matrix(con_ref)
