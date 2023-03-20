module Certificate

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import MultivariateBases as MB
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
function cone end
function gram_basis end
# FIXME currently, this returns `MB.MonomialBasis` instead of `MB.MonomialBasis{MT, MVT}`
function gram_basis_type end
function reduced_polynomial end
function ideal_certificate end
function preprocessed_domain end

# Only for PreorderIndex
function preorder_indices end
function generator end
function multiplier_basis end
function multiplier_basis_type end

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

cone(certificate::Putinar) = certificate.cone

struct WithVariables{S, V}
    inner::S
    variables::V
end

function MP.variables(v::WithVariables)
    return v.variables
end
function MP.monomials(v::WithVariables)
    return MP.monomials(v.inner)
end

_merge_sorted(a::Vector, ::Tuple{}) = a
function _merge_sorted(a::Vector, b::Vector)
    vars = sort!(vcat(a, b), rev = true)
    unique!(vars)
    return vars
end
_merge_sorted(a::Tuple{}, ::Tuple{}) = a
_merge_sorted(a::Tuple, ::Tuple{}) = a
_merge_sorted(::Tuple{}, b::Tuple) = b
function _merge_sorted(a::Tuple, b::Tuple)
    v = first(a)
    w = first(b)
    if v == w
        return (v, _merge_sorted(Base.tail(a), Base.tail(b))...)
    elseif v > w
        return (v, _merge_sorted(Base.tail(a), b)...)
    else
        return (w, _merge_sorted(a, Base.tail(b))...)
    end
end

_vars(::SemialgebraicSets.FullSpace) = tuple()
_vars(x) = MP.variables(x)

function with_variables(inner, outer)
    return WithVariables(inner, _merge_sorted(_vars(inner), _vars(outer)))
end

function preprocessed_domain(::Putinar, domain::BasicSemialgebraicSet, p)
    return with_variables(domain, p)
end

function preorder_indices(::Putinar, domain::WithVariables)
    return map(PreorderIndex, eachindex(domain.inner.p))
end

function maxdegree_gram_basis(B::Type, variables, maxdegree::Int)
    return MB.maxdegree_basis(B, variables, div(maxdegree, 2))
end
multiplier_maxdegree(maxdegree, q) = maxdegree - MP.maxdegree(q)
function multiplier_basis(certificate::Putinar, index::PreorderIndex, domain::WithVariables)
    q = domain.inner.p[index.value]
    vars = sort!([domain.variables..., MP.variables(q)...], rev = true)
    unique!(vars)
    return maxdegree_gram_basis(certificate.basis, vars, multiplier_maxdegree(certificate.maxdegree, q))
end
function multiplier_basis_type(::Type{Putinar{IC, CT, BT}}) where {IC, CT, BT}
    return BT
end

function generator(::Putinar, index::PreorderIndex, domain::WithVariables)
    return domain.inner.p[index.value]
end

ideal_certificate(certificate::Putinar) = certificate.ideal_certificate
ideal_certificate(::Type{<:Putinar{IC}}) where {IC} = IC

SumOfSquares.matrix_cone_type(::Type{<:Putinar{IC, CT}}) where {IC, CT} = SumOfSquares.matrix_cone_type(CT)

######################
# Ideal certificates #
######################

abstract type SimpleIdealCertificate{CT, BT} <: AbstractIdealCertificate end
reduced_polynomial(::SimpleIdealCertificate, poly, domain) = poly

cone(certificate::SimpleIdealCertificate) = certificate.cone
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
function gram_basis(certificate::MaxDegree, poly)
    return maxdegree_gram_basis(certificate.basis, MP.variables(poly), certificate.maxdegree)
end
function gram_basis_type(::Type{MaxDegree{CT, BT}}) where {CT, BT}
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
function gram_basis(certificate::FixedBasis, poly)
    return certificate.basis
end
function gram_basis_type(::Type{FixedBasis{CT, BT}}) where {CT, BT}
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
function gram_basis(certificate::Newton{CT, B}, poly) where {CT, B}
    return MB.basis_covering_monomials(B, monomials_half_newton_polytope(MP.monomials(poly), certificate.variable_groups))
end
function gram_basis_type(::Type{<:Newton{CT, BT}}) where {CT, BT}
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

function reduced_polynomial(::Remainder, poly, domain)
    return convert(typeof(poly), rem(poly, ideal(domain)))
end

function gram_basis(certificate::Remainder, poly)
    return gram_basis(certificate.gram_certificate, poly)
end
function gram_basis_type(::Type{Remainder{GCT}}) where GCT
    return gram_basis_type(GCT)
end

cone(certificate::Remainder) = cone(certificate.gram_certificate)
SumOfSquares.matrix_cone_type(::Type{Remainder{GCT}}) where {GCT} = SumOfSquares.matrix_cone_type(GCT)
zero_basis(certificate::Remainder) = zero_basis(certificate.gram_certificate)
zero_basis_type(::Type{Remainder{GCT}}) where {GCT} = zero_basis_type(GCT)

include("Sparsity/Sparsity.jl")
using .Sparsity: SignSymmetry, ChordalCompletion, ClusterCompletion
export Sparsity, SignSymmetry, ChordalCompletion, ClusterCompletion

include("Symmetry/Symmetry.jl")
export Symmetry

end
