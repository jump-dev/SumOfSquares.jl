export Sparsity, NoSparsity, VariableSparsity, MonomialSparsity, SignSymmetry
abstract type Sparsity end
struct NoSparsity <: Sparsity end

include("xor_space.jl")
include("sign.jl")
include("variable_sparsity.jl")
include("monomial_sparsity.jl")

struct ChordalPutinar{CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: AbstractPreorderCertificate
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end

get(certificate::ChordalPutinar, ::Cone) = certificate.cone

struct ChordalDomain{S, V}
    domain::S
    cliques::Vector{V}
end

function get(::ChordalPutinar, ::PreprocessedDomain, domain::BasicSemialgebraicSet, p)
    H, cliques = chordal_csp_graph(p, domain)
    return ChordalDomain(domain, cliques)
end

function get(::ChordalPutinar, index::PreorderIndices, domain::ChordalDomain)
    return map(PreorderIndex, eachindex(domain.domain.p))
end

function get(certificate::ChordalPutinar, ::MultiplierBasis, index::PreorderIndex, domain::ChordalDomain)
    q = domain.domain.p[index.value]
    maxdegree_s2 = certificate.maxdegree - MP.maxdegree(q)
    # If maxdegree_s2 is odd, `div(maxdegree_s2, 2)` would make s^2 have degree up to maxdegree_s2-1
    # for this reason, we take `div(maxdegree_s2 + 1, 2)` so that s^2 have degree up to maxdegree_s2+1
    return [maxdegree_gram_basis(certificate.basis, clique, maxdegree_s2 + 1) for clique in domain.cliques if variables(q) ⊆ clique]
end
function get(::Type{ChordalPutinar{CT, BT}}, ::MultiplierBasisType) where {CT, BT}
    return Vector{BT}
end

function get(::ChordalPutinar, ::Generator, index::PreorderIndex, domain::ChordalDomain)
    return domain.domain.p[index.value]
end

get(certificate::ChordalPutinar, ::IdealCertificate) = ChordalIdeal(MonomialSparsity(), certificate.cone, certificate.basis, certificate.maxdegree)
get(::Type{<:ChordalPutinar{CT, BT}}, ::IdealCertificate) where {CT, BT} = ChordalIdeal{CT, BT}

SumOfSquares.matrix_cone_type(::Type{<:ChordalPutinar{CT}}) where {CT} = SumOfSquares.matrix_cone_type(CT)

struct ChordalIdeal{S <: Sparsity, CT <: SumOfSquares.SOSLikeCone, BT <: MB.AbstractPolynomialBasis} <: SimpleIdealCertificate{CT, BT}
    sparsity::S
    cone::CT
    basis::Type{BT}
    maxdegree::Int
end
function sparsity(poly::MP.AbstractPolynomial, sp::VariableSparsity, basis, maxdegree)
    H, cliques = chordal_csp_graph(poly, FullSpace())
    return map(cliques) do clique
        return maxdegree_gram_basis(basis, clique, maxdegree)
    end
end
function monomial_sparsity()
end
function sparsity(poly::MP.AbstractPolynomial, sp::Union{SignSymmetry, MonomialSparsity}, basis::Type{<:MB.MonomialBasis}=MB.MonomialBasis, maxdegree=nothing)
    return MB.MonomialBasis.(sparsity(monomials(poly), sp))
end
function get(certificate::ChordalIdeal, ::GramBasis, poly)
    return sparsity(poly, certificate.sparsity, certificate.basis, certificate.maxdegree)
end
function get(::Type{ChordalIdeal{S, CT, BT}}, ::GramBasisType) where {S, CT, BT}
    return Vector{BT}
end
