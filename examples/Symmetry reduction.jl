# # Lyapunov Function Search

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Lyapunov Function Search.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Lyapunov Function Search.ipynb)
# **Adapted from**: SOSTOOLS' SOSDEMO2 (See Section 4.2 of [SOSTOOLS User's Manual](http://sysos.eng.ox.ac.uk/sostools/sostools.pdf))

using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl#tw/ex_sos"

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

function SymbolicWedderburn.OnPoints(basis::Union{MonomialVector, AbstractVector{<:Monomial}})
    basis_exps = Vector{Vector{Int}}(undef, length(basis))
    basis_dict = Dict{Vector{Int}, Int}()
    sizehint!(basis_dict, length(basis))

    for (i, b) in enumerate(basis)
        e = MP.exponents(b) # so that we allocate exponents only once
        basis_exps[i] = e
        basis_dict[e] = i
    end

    return SymbolicWedderburn.OnPoints(basis_exps, basis_dict)
end
function symmetry_adapted_basis(G, mvec)
    chars_vars = SymbolicWedderburn.characters_dixon(G)

    chars = chars_vars
    basis = mvec

    @assert all(χ.inv_of == first(chars).inv_of for χ in chars)

    induced_action = SymbolicWedderburn.OnPoints(basis)
    ccG_large = induced_action.(conjugacy_classes(first(chars)))

    let ccls = ccG_large, large_gens = induced_action.(gens(G)) # double check:
        G_large = PermGroup(large_gens)
        ccG_large = conjugacy_classes(G_large)
        @assert all(Set.(collect.(ccG_large)) .== Set.(collect.(ccls)))
    end

    chars_mvec = [SymbolicWedderburn.Character(values(χ), χ.inv_of, ccG_large) for χ in chars]

    vr_chars = SymbolicWedderburn.real_vchars(chars_mvec)
    U = filter!(x->all(!iszero, x), [SymbolicWedderburn.matrix_projection(χ) for χ in vr_chars])
    R = map(U) do c_u
        u = last(c_u)
        if all(isreal, u)
            image_coeffs, pivots = SymbolicWedderburn.row_echelon_form(float.(u))
            dim = length(pivots)
            image_coeffs[1:dim, :]
        else
            throw("Not Implemented")
        end
    end

    return map(R) do Ri
        FixedPolynomialBasis(Ri * mvec)
    end
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
    return symmetry_adapted_basis(cert.group, basis.monomials)
end
G = PermGroup([perm"(1,2,3,4)"])

# We can use this certificate as follows:

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, sum(x) + sum(x.^2) - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G))
optimize!(model)
@test value(t) ≈ -1 #src
value(t)

# We indeed find `-1`, let's verify that symmetry was exploited:

@test length(gram_matrix(con_ref).sub_gram_matrices) == 3 #src
gram_matrix(con_ref)
