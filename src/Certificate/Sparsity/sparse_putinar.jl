export Sparsity, NoSparsity, VariableSparsity, MonomialSparsity, SignSymmetry
abstract type Sparsity end
struct NoSparsity <: Sparsity end

include("xor_space.jl")
include("sign.jl")
include("variable_sparsity.jl")
include("monomial_sparsity.jl")

"""
    struct SparsePreorder{S <: Sparsity, C <: AbstractPreorderCertificate} <: AbstractPreorderCertificate
        sparsity::S
        certificate::C
    end

Same certificate as `C` except that the Sum-of-Squares polynomials `σ_i`
are modelled as a sum of Sum-of-Squares polynomials with smaller basis
using the sparsity reduction `sparsity`.
"""
struct SparsePreorder{S <: Sparsity, C <: SumOfSquares.Certificate.AbstractPreorderCertificate} <: SumOfSquares.Certificate.AbstractPreorderCertificate
    sparsity::S
    certificate::C
end

SumOfSquares.Certificate.get(certificate::SparsePreorder, attr::SumOfSquares.Certificate.Cone) = get(certificate.certificate, attr)

struct SparseDomain{S, P, B}
    domain::S
    processed::P
    bases::Vector{Vector{B}}
end

function SumOfSquares.Certificate.get(certificate::SparsePreorder, attr::SumOfSquares.Certificate.PreprocessedDomain, domain::BasicSemialgebraicSet, p)
    basis, preorder_bases = sparsity(p, domain, certificate.sparsity, certificate.certificate)
    return SparseDomain(domain, get(certificate.certificate, attr, domain, p), preorder_bases)
end

function SumOfSquares.Certificate.get(certificate::SparsePreorder, attr::SumOfSquares.Certificate.PreorderIndices, domain::SparseDomain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, domain.processed)
end

function SumOfSquares.Certificate.get(::SparsePreorder, ::SumOfSquares.Certificate.MultiplierBasis, index::SumOfSquares.Certificate.PreorderIndex, domain::SparseDomain)
    return domain.bases[index.value]
end
function SumOfSquares.Certificate.get(::Type{SparsePreorder{S, C}}, attr::SumOfSquares.Certificate.MultiplierBasisType) where {S, C}
    return Vector{SumOfSquares.Certificate.get(C, attr)}
end

function SumOfSquares.Certificate.get(certificate::SparsePreorder, attr::SumOfSquares.Certificate.Generator, index::SumOfSquares.Certificate.PreorderIndex, domain::SparseDomain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, index, domain.processed)
end

SumOfSquares.Certificate.get(certificate::SparsePreorder, attr::SumOfSquares.Certificate.IdealCertificate) = SparseIdeal(certificate.sparsity, SumOfSquares.Certificate.get(certificate.certificate, attr))
SumOfSquares.Certificate.get(::Type{<:SparsePreorder{S, C}}, attr::SumOfSquares.Certificate.IdealCertificate) where {S, C} = SparseIdeal{S, SumOfSquares.Certificate.get(C, attr)}

SumOfSquares.matrix_cone_type(::Type{SparsePreorder{S, C}}) where {S, C} = SumOfSquares.matrix_cone_type(C)

"""
    struct SparseIdeal{S <: Sparsity, C <: AbstractIdealCertificate} <: AbstractIdealCertificate
        sparsity::S
        certificate::C
    end

Same certificate as `C` except that the Sum-of-Squares polynomial `σ`
is modelled as a sum of Sum-of-Squares polynomials with smaller basis
using the sparsity reduction `sparsity`.
"""
struct SparseIdeal{S <: Sparsity, C <: SumOfSquares.Certificate.AbstractIdealCertificate} <: SumOfSquares.Certificate.AbstractIdealCertificate
    sparsity::S
    certificate::C
end

function SparseIdeal(sp::VariableSparsity, basis, cone, maxdegree::Nothing, newton_polytope)
    error("`maxdegree` cannot be `nothing` when `sparsity` is `VariableSparsity`.")
end
function SparseIdeal(sp::VariableSparsity, basis, cone, maxdegree::Integer, newton_polytope)
    return SparseIdeal(sp, SumOfSquares.Certificate.MaxDegree(cone, basis, maxdegree))
end
function SparseIdeal(sp::Union{MonomialSparsity, SignSymmetry}, basis, cone, maxdegree, newton_polytope)
    return SparseIdeal(sp, SumOfSquares.Certificate.Newton(cone, basis, newton_polytope))
end

function sparsity(poly::MP.AbstractPolynomial, ::VariableSparsity, certificate::SumOfSquares.Certificate.MaxDegree)
    H, cliques = chordal_csp_graph(poly, FullSpace())
    return map(cliques) do clique
        return SumOfSquares.Certificate.maxdegree_gram_basis(certificate.basis, clique, certificate.maxdegree)
    end
end
function sparsity(monos, sp::Union{SignSymmetry, MonomialSparsity}, gram_basis::MB.MonomialBasis)
    return MB.MonomialBasis.(sparsity(monos, sp, gram_basis.monomials))
end
function sparsity(poly::MP.AbstractPolynomial, sp::Union{SignSymmetry, MonomialSparsity}, certificate::SumOfSquares.Certificate.AbstractIdealCertificate)
    return sparsity(monomials(poly), sp, get(certificate, GramBasis(), poly))
end
function SumOfSquares.Certificate.get(certificate::SparseIdeal, ::SumOfSquares.Certificate.GramBasis, poly)
    return sparsity(poly, certificate.sparsity, certificate.certificate)
end
function SumOfSquares.Certificate.get(::Type{SparseIdeal{S, C}}, attr::SumOfSquares.Certificate.GramBasisType) where {S, C}
    return Vector{<:SumOfSquares.Certificate.get(C, attr)}
end
function SumOfSquares.Certificate.get(certificate::SparseIdeal, attr::SumOfSquares.Certificate.ReducedPolynomial, poly, domain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, poly, domain)
end
function SumOfSquares.Certificate.get(certificate::SparseIdeal, attr::SumOfSquares.Certificate.Cone)
    return SumOfSquares.Certificate.get(certificate.certificate, attr)
end
SumOfSquares.matrix_cone_type(::Type{SparseIdeal{S, C}}) where {S, C} = SumOfSquares.matrix_cone_type(C)

SumOfSquares.Certificate.zero_basis(certificate::SparseIdeal) = SumOfSquares.Certificate.zero_basis(certificate.certificate)
SumOfSquares.Certificate.zero_basis_type(::Type{SparseIdeal{S, C}}) where {S, C} = SumOfSquares.Certificate.zero_basis_type(C)
