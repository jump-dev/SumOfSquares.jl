export Sparsity, NoSparsity, VariableSparsity, MonomialSparsity, SignSymmetry
abstract type Sparsity end
struct NoSparsity <: Sparsity end

include("xor_space.jl")
include("sign.jl")
include("variable_sparsity.jl")
include("monomial_sparsity.jl")

struct SparsePreorder{S <: Sparsity, C <: AbstractPreorderCertificate} <: AbstractPreorderCertificate
    sparsity::S
    certificate::C
end

get(certificate::SparsePreorder, attr::Cone) = get(certificate.certificate, attr)

struct SparseDomain{S, P, B}
    domain::S
    processed::P
    bases::Vector{Vector{B}}
end

function get(certificate::SparsePreorder, attr::PreprocessedDomain, domain::BasicSemialgebraicSet, p)
    basis, preorder_bases = sparsity(p, domain, certificate.sparsity, certificate.certificate)
    return SparseDomain(domain, get(certificate.certificate, attr, domain, p), preorder_bases)
end

function get(certificate::SparsePreorder, attr::PreorderIndices, domain::SparseDomain)
    return get(certificate.certificate, attr, domain.processed)
end

function get(::SparsePreorder, ::MultiplierBasis, index::PreorderIndex, domain::SparseDomain)
    return domain.bases[index.value]
end
function get(::Type{SparsePreorder{S, C}}, attr::MultiplierBasisType) where {S, C}
    return Vector{get(C, attr)}
end

function get(certificate::SparsePreorder, attr::Generator, index::PreorderIndex, domain::SparseDomain)
    return get(certificate.certificate, attr, index, domain.processed)
end

get(certificate::SparsePreorder, attr::IdealCertificate) = SparseIdeal(certificate.sparsity, get(certificate.certificate, attr))
get(::Type{<:SparsePreorder{S, C}}, attr::IdealCertificate) where {S, C} = SparseIdeal{S, get(C, attr)}

SumOfSquares.matrix_cone_type(::Type{SparsePreorder{S, C}}) where {S, C} = SumOfSquares.matrix_cone_type(C)

struct SparseIdeal{S <: Sparsity, C <: AbstractIdealCertificate} <: AbstractIdealCertificate
    sparsity::S
    certificate::C
end

function SparseIdeal(sp::VariableSparsity, cone, basis, maxdegree::Nothing, newton_polytope)
    error("`maxdegree` cannot be `nothing` when `sparsity` is `VariableSparsity`.")
end
function SparseIdeal(sp::VariableSparsity, cone, basis, maxdegree::Integer, newton_polytope)
    return SparseIdeal(sp, MaxDegree(cone, basis, maxdegree))
end
function SparseIdeal(sp::Union{MonomialSparsity, SignSymmetry}, cone, basis, maxdegree, newton_polytope)
    return SparseIdeal(sp, Newton(cone, basis, newton_polytope))
end

function sparsity(poly::MP.AbstractPolynomial, ::VariableSparsity, certificate::MaxDegree)
    H, cliques = chordal_csp_graph(poly, FullSpace())
    return map(cliques) do clique
        return maxdegree_gram_basis(certificate.basis, clique, certificate.maxdegree)
    end
end
function sparsity(monos, sp::Union{SignSymmetry, MonomialSparsity}, gram_basis::MB.MonomialBasis)
    return MB.MonomialBasis.(sparsity(monos, sp, gram_basis.monomials))
end
function sparsity(poly::MP.AbstractPolynomial, sp::Union{SignSymmetry, MonomialSparsity}, certificate::AbstractIdealCertificate)
    return sparsity(monomials(poly), sp, get(certificate, GramBasis(), poly))
end
function get(certificate::SparseIdeal, ::GramBasis, poly)
    return sparsity(poly, certificate.sparsity, certificate.certificate)
end
function get(::Type{SparseIdeal{S, C}}, attr::GramBasisType) where {S, C}
    return Vector{<:get(C, attr)}
end
function get(certificate::SparseIdeal, attr::ReducedPolynomial, poly, domain)
    return get(certificate.certificate, attr, poly, domain)
end
function get(certificate::SparseIdeal, attr::Cone)
    return get(certificate.certificate, attr)
end
SumOfSquares.matrix_cone_type(::Type{SparseIdeal{S, C}}) where {S, C} = SumOfSquares.matrix_cone_type(C)

zero_basis(certificate::SparseIdeal) = zero_basis(certificate.certificate)
zero_basis_type(::Type{SparseIdeal{S, C}}) where {S, C} = zero_basis_type(C)
