using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl#enh/symmetry_adapted_basis"

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using DynamicPolynomials
@polyvar x[1:4]

poly = sum(x) + sum(x.^2)

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

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, sum(x) + sum(x.^2) - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G))
optimize!(model)
value(t)

gram_matrix(con_ref)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

