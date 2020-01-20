module Certificate

import MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
using SemialgebraicSets

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
struct GramBasisType <: Attribute end
struct ReducedPolynomial <: Attribute end
struct IdealCertificate <: Attribute end
struct PreprocessedDomain <: Attribute end

# Only for PreorderIndex
struct PreorderIndices <: Attribute end
struct Generator <: Attribute end
struct MultiplierBasis <: Attribute end
struct MultiplierBasisType <: Attribute end

# For set
#struct Monomials <: Attribute end


# PreorderCertificate
# IdealCertificate
# FullSpaceCertificate

abstract type AbstractCertificate end

abstract type AbstractPreorderCertificate <: AbstractCertificate end
abstract type AbstractIdealCertificate <: AbstractCertificate end

struct Putinar{IC <: AbstractIdealCertificate, CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: AbstractPreorderCertificate
    ideal_certificate::IC
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end

get(certificate::Putinar, ::Cone) = certificate.cone

struct DomainWithVariables{S, V}
    domain::S
    variables::V
end
function get(::Putinar, ::PreprocessedDomain, domain::BasicSemialgebraicSet, p)
    return DomainWithVariables(domain, MP.variables(p))
end

function get(::Putinar, index::PreorderIndices, domain::DomainWithVariables)
    return map(PreorderIndex, eachindex(domain.domain.p))
end

function get(certificate::Putinar, ::MultiplierBasis, index::PreorderIndex, domain::DomainWithVariables)
    q = domain.domain.p[index.value]
    vars = sort!([domain.variables..., MP.variables(q)...], rev = true)
    unique!(vars)
    maxdegree_s2 = certificate.maxdegree - MP.maxdegree(q)
    # If maxdegree_s2 is odd, `div(maxdegree_s2, 2)` would make s^2 have degree up to maxdegree_s2-1
    # for this reason, we take `div(maxdegree_s2 + 1, 2)` so that s^2 have degree up to maxdegree_s2+1
    return maxdegree_basis(certificate.basis, vars, maxdegree_s2 + 1)
end
function get(::Type{Putinar{IC, CT, BT}}, ::MultiplierBasisType) where {IC, CT, BT}
    return BT
end

function get(::Putinar, ::Generator, index::PreorderIndex, domain::DomainWithVariables)
    return domain.domain.p[index.value]
end

get(certificate::Putinar, ::IdealCertificate) = certificate.ideal_certificate
get(::Type{<:Putinar{IC}}, ::IdealCertificate) where {IC} = IC

SumOfSquares.matrix_cone_type(::Type{<:Putinar{IC, CT}}) where {IC, CT} = SumOfSquares.matrix_cone_type(CT)

######################
# Ideal certificates #
######################

abstract type SimpleIdealCertificate{CT, BT} <: AbstractIdealCertificate end
get(::SimpleIdealCertificate, ::ReducedPolynomial, poly, domain) = poly

get(certificate::SimpleIdealCertificate, ::Cone) = certificate.cone
SumOfSquares.matrix_cone_type(::Type{<:SimpleIdealCertificate{CT}}) where {CT} = SumOfSquares.matrix_cone_type(CT)
zero_basis(certificate::SimpleIdealCertificate) = certificate.basis
zero_basis_type(::Type{<:SimpleIdealCertificate{CT, BT}}) where {CT, BT} = BT

struct MaxDegree{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end
function maxdegree_basis(::Type{MB.MonomialBasis}, variables, maxdegree::Int)
    return MB.MonomialBasis(MP.monomials(variables, 0:div(maxdegree, 2)))
end
function get(certificate::MaxDegree, ::GramBasis, poly) where CT
    return maxdegree_basis(certificate.basis, MP.variables(poly), certificate.maxdegree)
end
function get(::Type{MaxDegree{CT, BT}}, ::GramBasisType) where {CT, BT}
    return BT
end

struct Newton{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis, NPT <: Tuple} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::Type{BT}
    variable_groups::NPT
end
function get(certificate::Newton{CT, MB.MonomialBasis}, ::GramBasis, poly) where CT
    return MB.MonomialBasis(monomials_half_newton_polytope(MP.monomials(poly), certificate.variable_groups))
end
function get(::Type{<:Newton{CT, BT}}, ::GramBasisType) where {CT, BT}
    return BT
end

struct Remainder{GCT<:AbstractIdealCertificate} <: AbstractIdealCertificate
    gram_certificate::GCT
end

function get(::Remainder, ::ReducedPolynomial, poly, domain)
    return convert(typeof(poly), rem(poly, ideal(domain)))
end

function get(certificate::Remainder, attr::GramBasis, poly)
    return get(certificate.gram_certificate, attr, poly)
end
function get(::Type{Remainder{GCT}}, attr::GramBasisType) where GCT
    return get(GCT, attr)
end

get(certificate::Remainder, attr::Cone) = get(certificate.gram_certificate, attr)
SumOfSquares.matrix_cone_type(::Type{Remainder{GCT}}) where {GCT} = SumOfSquares.matrix_cone_type(GCT)
zero_basis(certificate::Remainder) = zero_basis(certificate.gram_certificate)
zero_basis_type(::Type{Remainder{GCT}}) where {GCT} = zero_basis_type(GCT)

include("csp/ChordalExtensionGraph.jl")
include("csp/sparse_putinar.jl")

end
