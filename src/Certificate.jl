module Certificate

import MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets

# TODO replace by MultivariateBases
using PolyJuMP

using SumOfSquares

include("newton_polytope.jl")

# Certificate              denominator * p = ...                                      domain
# Polya       (1 + x1 + ... + xn)^(2d) * p = ∑ λ_α x.^α                  λ_α ∈ R_+    R_+^n
# Putinar                                p = σ_0 + ∑ σ_i g_i             σ_i SOS      g_i ≥ 0
# Krivine                                p = ∑ λ_{αβ} g.^α (1 .- g).^β   λ_{αβ} ≥ 0   0 ≤ g_i ≤ 1
# Schmüdgen                              p = ∑ σ_α g^α                   σ_α SOS      g_i ≥ 0
# Hilbert                            d * p = σ                           d, σ SOS     R^n

struct PreorderIndex
    value::Int
end
struct IdealIndex
    value::Int
end

const Index = Union{PreorderIndex, IdealIndex}

abstract type Attribute end

# For get
struct Cone <: Attribute end
struct GramBasis <: Attribute end
struct ReducedPolynomial <: Attribute end
struct IdealCertificate <: Attribute end

# Only for PreorderIndex
struct PreorderIndices <: Attribute end
struct Generator <: Attribute end
struct MultiplierBasis <: Attribute end

# For set
#struct Monomials <: Attribute end


# PreorderCertificate
# IdealCertificate
# FullSpaceCertificate

abstract type AbstractCertificate end

abstract type AbstractPreorderCertificate <: AbstractCertificate end
abstract type AbstractIdealCertificate <: AbstractCertificate end

struct Putinar{IC <: AbstractIdealCertificate, CT <: SumOfSquares.SOSLikeCone, BT <: PolyJuMP.AbstractPolynomialBasis} <: AbstractPreorderCertificate
    ideal_certificate::IC
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end

function get(::Putinar, index::PreorderIndices, domain::BasicSemialgebraicSet)
    return map(PreorderIndex, eachindex(domain.p))
end

function get(certificate::Putinar, ::MultiplierBasis, index::PreorderIndex, domain::BasicSemialgebraicSet, p)
    q = domain.p[index.value]
    maxdegree_s2 = certificate.maxdegree - MP.maxdegree(q)
    # If maxdegree_s2 is odd, `div(maxdegree_s2, 2)` would make s^2 have degree up to maxdegree_s2-1
    # for this reason, we take `div(maxdegree_s2 + 1, 2)` so that s^2 have degree up to maxdegree_s2+1
    maxdegree_s = div(maxdegree_s2 + 1, 2)
    vars = sort!([MP.variables(p)..., MP.variables(q)...], rev = true)
    unique!(vars)
    return MP.monomials(vars, 0:maxdegree_s)
end

function get(::Putinar, ::Generator, index::PreorderIndex, domain::BasicSemialgebraicSet)
    return domain.p[index.value]
end

get(certificate::Putinar, ::IdealCertificate) = certificate.ideal_certificate
get(::Type{<:Putinar{IC}}, ::IdealCertificate) where {IC} = IC

SumOfSquares.matrix_cone_type(::Type{<:Putinar{IC, CT}}) where {IC, CT} = SumOfSquares.matrix_cone_type(CT)

######################
# Ideal certificates #
######################

struct MaxDegree{CT <: SumOfSquares.SOSLikeCone, BT <: PolyJuMP.AbstractPolynomialBasis} <: AbstractIdealCertificate
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end
function get(certificate::MaxDegree, ::GramBasis, poly)
    return monomials(MP.variables(poly), 0:div(certificate.maxdegree, 2))
end

struct Newton{CT <: SumOfSquares.SOSLikeCone, BT <: PolyJuMP.AbstractPolynomialBasis, NPT <: Tuple} <: AbstractIdealCertificate
    cone::CT
    basis::Type{BT}
    variable_groups::NPT
end
function get(certificate::Newton, ::GramBasis, poly)
    return monomials_half_newton_polytope(MP.monomials(poly), certificate.variable_groups)
end

get(::Union{MaxDegree, Newton}, ::ReducedPolynomial, poly, domain) = poly

SumOfSquares.matrix_cone_type(::Type{<:Union{MaxDegree{CT}, Newton{CT}}}) where {CT} = SumOfSquares.matrix_cone_type(CT)
zero_basis(certificate::Union{MaxDegree, Newton}) = certificate.basis
zero_basis_type(::Type{<:Union{MaxDegree{CT, BT}, Newton{CT, BT}}}) where {CT, BT} = BT

struct Remainder{GCT<:AbstractIdealCertificate} <: AbstractIdealCertificate
    gram_certificate::GCT
end

function get(::Remainder, ::ReducedPolynomial, poly, domain)
    return convert(typeof(poly), rem(poly, ideal(domain)))
end

function get(certificate::Remainder, attr::GramBasis, poly)
    return get(certificate.gram_certificate, attr, poly)
end

SumOfSquares.matrix_cone_type(::Type{Remainder{GCT}}) where {GCT} = SumOfSquares.matrix_cone_type(GCT)
zero_basis(certificate::Remainder) = zero_basis(certificate.gram_certificate)
zero_basis_type(::Type{Remainder{GCT}}) where {GCT} = zero_basis_type(GCT)

end
