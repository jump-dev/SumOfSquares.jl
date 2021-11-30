module Certificate

import MutableArithmetics
const MA = MutableArithmetics
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
# FIXME currently, this returns `MB.MonomialBasis` instead of `MB.MonomialBasis{MT, MVT}`
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

"""
    struct Putinar{IC <: AbstractIdealCertificate, CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: AbstractPreorderCertificate
        ideal_certificate::IC
        cone::CT
        basis::Type{BT}
        maxdegree::Int
    end

The `Putinar` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`g_i(x) >= 0` and `h_i(x) = 0` by exhibiting Sum-of-Squares polynomials `σ_i(x)`
such that `p(x) - ∑ σ_i(x) g_i(x)` is guaranteed to be nonnegativity for all `x`
such that `h_i(x) = 0`.
The polynomials `σ_i(x)` are search over `cone` with a basis of type `basis` such that
the degree of `σ_i(x) g_i(x)` does not exceed `maxdegree`.
"""
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

function get(::Putinar, ::PreorderIndices, domain::DomainWithVariables)
    return map(PreorderIndex, eachindex(domain.domain.p))
end

function maxdegree_gram_basis(B::Type, variables, maxdegree::Int)
    return MB.maxdegree_basis(B, variables, div(maxdegree, 2))
end
multiplier_maxdegree(maxdegree, q) = maxdegree - MP.maxdegree(q)
function get(certificate::Putinar, ::MultiplierBasis, index::PreorderIndex, domain::DomainWithVariables)
    q = domain.domain.p[index.value]
    vars = sort!([domain.variables..., MP.variables(q)...], rev = true)
    unique!(vars)
    return maxdegree_gram_basis(certificate.basis, vars, multiplier_maxdegree(certificate.maxdegree, q))
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

# TODO return something else when `PolyJuMP` support other bases.
zero_basis(certificate::SimpleIdealCertificate) = MB.MonomialBasis
zero_basis_type(::Type{<:SimpleIdealCertificate{CT, BT}}) where {CT, BT} = MB.MonomialBasis

"""
    struct MaxDegree{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
        cone::CT
        basis::Type{BT}
        maxdegree::Int
    end

The `MaxDegree` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with a basis of type `basis` such that
the degree of `σ(x)` does not exceed `maxdegree`.
"""
struct MaxDegree{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end
function get(certificate::MaxDegree, ::GramBasis, poly) where CT
    return maxdegree_gram_basis(certificate.basis, MP.variables(poly), certificate.maxdegree)
end
function get(::Type{MaxDegree{CT, BT}}, ::GramBasisType) where {CT, BT}
    return BT
end

"""
    struct FixedBasis{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
        cone::CT
        basis::BT
    end

The `FixedBasis` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with basis `basis`.
"""
struct FixedBasis{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::BT
end
function get(certificate::FixedBasis, ::GramBasis, poly) where CT
    return certificate.basis
end
function get(::Type{FixedBasis{CT, BT}}, ::GramBasisType) where {CT, BT}
    return BT
end

"""
    struct Newton{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis, NPT <: Tuple} <: SimpleIdealCertificate{CT, BT}
        cone::CT
        basis::Type{BT}
        variable_groups::NPT
    end

The `Newton` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with a basis of type `basis`
chosen using the multipartite Newton polytope with parts `variable_groups`.
If `variable_groups = tuple()` then it falls back to the classical Newton polytope
with all variables in the same part.
"""
struct Newton{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis, NPT <: Tuple} <: SimpleIdealCertificate{CT, BT}
    cone::CT
    basis::Type{BT}
    variable_groups::NPT
end
function get(certificate::Newton{CT, B}, ::GramBasis, poly) where {CT, B}
    return MB.basis_covering_monomials(B, monomials_half_newton_polytope(MP.monomials(poly), certificate.variable_groups))
end
function get(::Type{<:Newton{CT, BT}}, ::GramBasisType) where {CT, BT}
    return BT
end

"""
    struct Remainder{GCT<:AbstractIdealCertificate} <: AbstractIdealCertificate
        gram_certificate::GCT
    end

The `Remainder` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by guaranteeing the remainder of `p(x)` modulo the ideal generated by
`⟨h_i⟩` to be nonnegative for all `x` such that `h_i(x) = 0` using the certificate
`gram_certificate`.
For instance, if `gram_certificate` is [`SumOfSquares.Certificate.Newton`](@ref),
then the certificate `Remainder(gram_certificate)` will take the remainder before
computing the Newton polytope hence might generate a much smaller Newton polytope
hence a smaller basis and smaller semidefinite program.
However, this then corresponds to a lower degree of the hierarchy which might
be insufficient to find a certificate.
"""
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

include("Sparsity/Sparsity.jl")
using .Sparsity: SignSymmetry, ChordalCompletion, ClusterCompletion
export Sparsity, SignSymmetry, ChordalCompletion, ClusterCompletion

include("Symmetry/Symmetry.jl")
export Symmetry

end
