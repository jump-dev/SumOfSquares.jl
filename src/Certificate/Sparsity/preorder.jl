"""
    struct Sparsity.Preorder{S <: Sparsity.Pattern, C <: SumOfSquares.Certificate.AbstractPreorderCertificate} <: SumOfSquares.Certificate.AbstractPreorderCertificate
        sparsity::S
        certificate::C
    end

Same certificate as `C` except that the Sum-of-Squares polynomials `σ_i`
are modelled as a sum of Sum-of-Squares polynomials with smaller basis
using the sparsity reduction `sparsity`.
"""
struct Preorder{S <: Pattern, C <: SumOfSquares.Certificate.AbstractPreorderCertificate} <: SumOfSquares.Certificate.AbstractPreorderCertificate
    sparsity::S
    certificate::C
end

SumOfSquares.Certificate.get(certificate::Preorder, attr::SumOfSquares.Certificate.Cone) = SumOfSquares.Certificate.get(certificate.certificate, attr)

struct Domain{S, P, B}
    domain::S
    processed::P
    bases::Vector{Vector{B}}
end

function SumOfSquares.Certificate.get(certificate::Preorder, attr::SumOfSquares.Certificate.PreprocessedDomain, domain::SemialgebraicSets.BasicSemialgebraicSet, p)
    basis, Preorder_bases = sparsity(p, domain, certificate.sparsity, certificate.certificate)
    return Domain(domain, SumOfSquares.Certificate.get(certificate.certificate, attr, domain, p), Preorder_bases)
end

function SumOfSquares.Certificate.get(certificate::Preorder, attr::SumOfSquares.Certificate.PreorderIndices, domain::Domain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, domain.processed)
end

function SumOfSquares.Certificate.get(::Preorder, ::SumOfSquares.Certificate.MultiplierBasis, index::SumOfSquares.Certificate.PreorderIndex, domain::Domain)
    return domain.bases[index.value]
end
function SumOfSquares.Certificate.get(::Type{Preorder{S, C}}, attr::SumOfSquares.Certificate.MultiplierBasisType) where {S, C}
    return Vector{SumOfSquares.Certificate.get(C, attr)}
end

function SumOfSquares.Certificate.get(certificate::Preorder, attr::SumOfSquares.Certificate.Generator, index::SumOfSquares.Certificate.PreorderIndex, domain::Domain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, index, domain.processed)
end

SumOfSquares.Certificate.get(certificate::Preorder, attr::SumOfSquares.Certificate.IdealCertificate) = Ideal(certificate.sparsity, SumOfSquares.Certificate.get(certificate.certificate, attr))
SumOfSquares.Certificate.get(::Type{<:Preorder{S, C}}, attr::SumOfSquares.Certificate.IdealCertificate) where {S, C} = Ideal{S, SumOfSquares.Certificate.get(C, attr)}

SumOfSquares.matrix_cone_type(::Type{Preorder{S, C}}) where {S, C} = SumOfSquares.matrix_cone_type(C)
