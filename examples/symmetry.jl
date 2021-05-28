using Pkg
pkg"add https://github.com/kalmarek/SymbolicWedderburn.jl#bl/sos"

import MutableArithmetics
const MA = MutableArithmetics
using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateBases
const MB = MultivariateBases

using GroupsCore
using SymbolicWedderburn
using PermutationGroups
using Cyclotomics
using SumOfSquares

function SymbolicWedderburn.decompose(k::MP.AbstractPolynomialLike, hom::SymbolicWedderburn.InducedActionHomomorphism)
    # correct only if features(hom) == monomials

    indcs = [hom[mono] for mono in MP.monomials(k)]
    coeffs = MP.coefficients(k)

    return indcs, coeffs
end

function SymbolicWedderburn.ExtensionHomomorphism(action::SymbolicWedderburn.Action, basis::MB.MonomialBasis)
    monos = collect(basis.monomials)
    mono_to_index = Dict(monos[i] => i for i in eachindex(monos))
    return SymbolicWedderburn.ExtensionHomomorphism(action, monos, mono_to_index)
end

function MP.polynomialtype(::Type{<:MB.AbstractPolynomialVectorBasis{PT}}, T::Type) where PT
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end

struct SymmetricIdeal{C, GT, AT} <: Certificate.AbstractIdealCertificate
    certificate::C
    group::GT
    action::AT
end
SumOfSquares.matrix_cone_type(::Type{<:SymmetricIdeal{C}}) where {C} = SumOfSquares.matrix_cone_type(C)
Certificate.get(::Type{<:SymmetricIdeal}, ::SumOfSquares.Certificate.GramBasisType) = Vector{Vector{MB.FixedPolynomialBasis}}
Certificate.zero_basis_type(::Type{<:SymmetricIdeal}) = MB.MonomialBasis
Certificate.zero_basis(::SymmetricIdeal) = MB.MonomialBasis
function Certificate.get(certificate::SymmetricIdeal, attr::Certificate.ReducedPolynomial, poly, domain)
    return Certificate.get(certificate.certificate, attr, poly, domain)
end
_type(::Type{MOI.PositiveSemidefiniteConeTriangle}) = Float64
_type(::Type{SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle}) = Complex{Float64}
function Certificate.get(cert::SymmetricIdeal, attr::Certificate.GramBasis, poly)
    basis = Certificate.get(cert.certificate, attr, poly)
    T = _type(SumOfSquares.matrix_cone_type(typeof(cert)))
    Rs, ms = SymbolicWedderburn.symmetry_adapted_basis(T, cert.group, basis, cert.action)
    return map(zip(Rs, ms)) do (R, m)
        F = convert(Matrix{T}, R)
        N = size(R, 1)
        d = div(N, m)
        if d > 1
            if m > 1
                ps = R * basis.monomials
                S = map(gens(cert.group)) do g
                    Si = Matrix{T}(undef, N, N)
                    for i in eachindex(ps)
                        p = ps[i]
                        q = SymbolicWedderburn.action(cert.action, g, p)
                        coefs = coefficients(q, basis.monomials)
                        col = row_echelon_linsolve(R, coefs)
                        Si[:, i] = col
                    end
                    return Si
                end
                U = ordered_block_diag(S, d)
            else
                U = Matrix(1.0I, N, N)
            end
            map(1:d) do i
                FixedPolynomialBasis((transpose(U[:, i:d:(i+d*(m-1))]) * F) * basis.monomials)
            end
        else
            [FixedPolynomialBasis(F * basis.monomials)]
        end
    end
end

include(joinpath(dirname(dirname(pathof(SumOfSquares))), "examples", "block_diag.jl"))
