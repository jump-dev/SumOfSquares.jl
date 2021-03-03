using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl#bl/nonperm"

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using SymbolicWedderburn
using PermutationGroups
using Cyclotomics
using SumOfSquares

function SymbolicWedderburn.ExtensionHomomorphism(basis::MB.MonomialBasis, action)
    monos = collect(basis.monomials)
    mono_to_index = Dict(monos[i] => i for i in eachindex(monos))
    return SymbolicWedderburn.ExtensionHomomorphism(monos, mono_to_index, action)
end

function MP.polynomialtype(::Type{<:MB.AbstractPolynomialVectorBasis{PT}}, T::Type) where PT
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end

struct SymmetricIdeal{CT, GT, AT} <: Certificate.AbstractIdealCertificate
    cone::CT
    group::GT
    action::AT
end
SumOfSquares.matrix_cone_type(::Type{<:SymmetricIdeal{CT}}) where {CT} = SumOfSquares.matrix_cone_type(CT)
Certificate.get(::Type{<:SymmetricIdeal}, ::SumOfSquares.Certificate.GramBasisType) = Vector{MB.FixedPolynomialBasis}
Certificate.zero_basis_type(::Type{<:SymmetricIdeal}) = MB.MonomialBasis
Certificate.zero_basis(::SymmetricIdeal) = MB.MonomialBasis
Certificate.get(::SymmetricIdeal, ::Certificate.ReducedPolynomial, poly, domain) = poly
function Certificate.get(cert::SymmetricIdeal, ::Certificate.GramBasis, poly)
    basis = Certificate.maxdegree_gram_basis(MB.MonomialBasis, MP.variables(poly), MP.maxdegree(poly))
    R = SymbolicWedderburn.symmetry_adapted_basis(Float64, cert.group, basis, cert.action)
    return map(R) do Ri
        FixedPolynomialBasis(convert(Matrix{Float64}, Ri) * basis.monomials)
    end
end
