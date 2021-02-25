# # Even reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Even reduction.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Even reduction.ipynb)

using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl#bl/nonperm"

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using Test #src
using DynamicPolynomials
@polyvar x

# We would like to find the minimum value of the polynomial

poly = x^4 - 2x^2

using SymbolicWedderburn
using PermutationGroups
using Cyclotomics
using SumOfSquares

# This is still a work in progress in SumOfSquares, so we need to define things here that will be moved inside SumOfSquares.jl once SymbolicWedderburn.jl is released:

struct ScaledPerm2{T, I} <: AbstractPerm
    indices::Vector{Pair{I, T}}
end
Base.one(p::ScaledPerm2{T}) where {T} = ScaledPerm2([i => one(T) for i in eachindex(p.indices)])
Base.:(==)(p::ScaledPerm2, q::ScaledPerm2) = p.indices == q.indices
Base.hash(p::ScaledPerm2, u::UInt64) = hash(p.indices, u)
SymbolicWedderburn.degree(p::ScaledPerm2) = length(p.indices)
Base.:^(i::Integer, p::ScaledPerm2) = p.indices[i]
function SymbolicWedderburn.add_inverse_permutation!(result, val, i::Int, j::Pair)
    result[i, j.first] += val / j.second
end

function permutation(ehom::SymbolicWedderburn.ExtensionHomomorphism{T}, els::Vector{T}) where {T}
    return Perm([ehom[el] for el in els])
end
function permutation(ehom::SymbolicWedderburn.ExtensionHomomorphism{<:AbstractMonomial}, terms::Vector{<:AbstractTerm})
    return ScaledPerm2([ehom[monomial(term)] => coefficient(term) for term in terms])
end
function (ehom::SymbolicWedderburn.ExtensionHomomorphism)(action)
    return permutation(ehom, [f^action for f in ehom.features])
end

function SymbolicWedderburn.ExtensionHomomorphism(basis::MB.MonomialBasis)
    monos = collect(basis.monomials)
    mono_to_index = Dict(monos[i] => i for i in eachindex(monos))
    return SymbolicWedderburn.ExtensionHomomorphism(monos, mono_to_index)
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
    R = SymbolicWedderburn.symmetry_adapted_basis(Float64, cert.group, basis)
    return map(R) do Ri
        FixedPolynomialBasis(convert(Matrix{Float64}, Ri) * basis.monomials)
    end
end


# We define the custom group as follows:

struct EvenOddAction <: GroupElem
    identity::Bool
end
PermutationGroups.order(a::EvenOddAction) = a.identity ? 1 : 2
Base.one(::EvenOddAction) = EvenOddAction(true)
Base.inv(a::EvenOddAction) = a
PermutationGroups.mul!(::EvenOddAction, a::EvenOddAction, b::EvenOddAction) = EvenOddAction(xor(a.identity, b.identity))
function Base.:^(mono::AbstractMonomial, a::EvenOddAction)
    if a.identity || iseven(MP.degree(mono))
        return 1 * mono
    else
        return -1 * mono
    end
end

struct EvenOddSymmetry <: Group
end
_orbit(cc) = PermutationGroups.Orbit(cc, Dict(a => nothing for a in cc))
SymbolicWedderburn.conjugacy_classes_orbit(::EvenOddSymmetry) = [_orbit([EvenOddAction(true)]), _orbit([EvenOddAction(false)])]

G = EvenOddSymmetry()

# We can exploit the symmetry as follows:

import CSDP
solver = CSDP.Optimizer
model = Model(solver)
@variable(model, t)
@objective(model, Max, t)
con_ref = @constraint(model, poly - t in SOSCone(), ideal_certificate = SymmetricIdeal(SOSCone(), G))
optimize!(model)
@show value(t)
value(t)

# We indeed find `-1`, let's verify that symmetry was exploited:

@test length(gram_matrix(con_ref).sub_gram_matrices) == 2 #src
@test gram_matrix(con_ref).sub_gram_matrices[1].basis.polynomials == [x] #src
@test gram_matrix(con_ref).sub_gram_matrices[2].basis.polynomials == [x^2, 1] #src
gram_matrix(con_ref)
