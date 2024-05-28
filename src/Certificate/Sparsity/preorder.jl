"""
    struct Sparsity.Preorder{S <: Sparsity.Pattern, C <: SumOfSquares.Certificate.AbstractPreorderCertificate} <: SumOfSquares.Certificate.AbstractPreorderCertificate
        sparsity::S
        certificate::C
    end

Same certificate as `C` except that the Sum-of-Squares polynomials `σ_i`
are modelled as a sum of Sum-of-Squares polynomials with smaller basis
using the sparsity reduction `sparsity`.
"""
struct Preorder{
    S<:Pattern,
    C<:SumOfSquares.Certificate.AbstractPreorderCertificate,
} <: SumOfSquares.Certificate.AbstractPreorderCertificate
    sparsity::S
    certificate::C
end

function SumOfSquares.Certificate.cone(certificate::Preorder)
    return SumOfSquares.Certificate.cone(certificate.certificate)
end

struct Domain{S,P,B}
    domain::S
    processed::P
    bases::Vector{Vector{B}}
end

function SumOfSquares.Certificate.preprocessed_domain(
    certificate::Preorder,
    domain::SemialgebraicSets.BasicSemialgebraicSet,
    p,
)
    basis, Preorder_bases =
        sparsity(MB.SubBasis{MB.Monomial}(MP.monomials(p)), domain, certificate.sparsity, certificate.certificate)
    return Domain(
        domain,
        SumOfSquares.Certificate.preprocessed_domain(
            certificate.certificate,
            domain,
            p,
        ),
        Preorder_bases,
    )
end

function SumOfSquares.Certificate.preorder_indices(
    certificate::Preorder,
    domain::Domain,
)
    return SumOfSquares.Certificate.preorder_indices(
        certificate.certificate,
        domain.processed,
    )
end

function SumOfSquares.Certificate.multiplier_basis(
    ::Preorder,
    index::SumOfSquares.Certificate.PreorderIndex,
    domain::Domain,
)
    return domain.bases[index.value]
end
function SumOfSquares.Certificate.multiplier_basis_type(
    ::Type{Preorder{S,C}},
    ::Type{M},
) where {S,C,M}
    return SumOfSquares.Certificate.multiplier_basis_type(C, M)
end

function SumOfSquares.Certificate.generator(
    certificate::Preorder,
    index::SumOfSquares.Certificate.PreorderIndex,
    domain::Domain,
)
    return SumOfSquares.Certificate.generator(
        certificate.certificate,
        index,
        domain.processed,
    )
end

function SumOfSquares.Certificate.ideal_certificate(certificate::Preorder)
    return Ideal(
        certificate.sparsity,
        SumOfSquares.Certificate.ideal_certificate(certificate.certificate),
    )
end
function SumOfSquares.Certificate.ideal_certificate(
    ::Type{<:Preorder{S,C}},
) where {S,C}
    return Ideal{S,SumOfSquares.Certificate.ideal_certificate(C)}
end

function SumOfSquares.matrix_cone_type(::Type{Preorder{S,C}}) where {S,C}
    return SumOfSquares.matrix_cone_type(C)
end
