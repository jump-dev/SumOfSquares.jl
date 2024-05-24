"""
    struct Sparsity.Ideal{S <: Sparsity.Pattern, C <: AbstractIdealCertificate} <: SumOfSquares.Certificate.AbstractIdealCertificate
        sparsity::S
        certificate::C
    end

Same certificate as `certificate` except that the Sum-of-Squares polynomial `σ`
is modelled as a sum of Sum-of-Squares polynomials with smaller bases
using the sparsity reduction `sparsity`.
"""
struct Ideal{S<:Pattern,C<:SumOfSquares.Certificate.AbstractIdealCertificate} <:
       SumOfSquares.Certificate.AbstractIdealCertificate
    sparsity::S
    certificate::C
end

function Ideal(sp::Variable, basis, cone, maxdegree::Nothing, newton_polytope)
    return error(
        "`maxdegree` cannot be `nothing` when `sparsity` is `Sparsity.Variable`.",
    )
end
function Ideal(sp::Variable, basis, cone, maxdegree::Integer, newton_polytope)
    return Ideal(sp, SumOfSquares.Certificate.MaxDegree(cone, basis, maxdegree))
end
function Ideal(
    sp::Union{Monomial,SignSymmetry},
    basis,
    cone,
    maxdegree,
    newton_polytope,
)
    return Ideal(
        sp,
        SumOfSquares.Certificate.Newton(cone, basis, newton_polytope),
    )
end

function sparsity(
    poly::MP.AbstractPolynomial,
    ::Variable,
    certificate::SumOfSquares.Certificate.MaxDegree,
)
    H, cliques = chordal_csp_graph(poly, SemialgebraicSets.FullSpace())
    return map(cliques) do clique
        return SumOfSquares.Certificate.maxdegree_gram_basis(
            certificate.basis,
            clique,
            certificate.maxdegree,
        )
    end
end
function sparsity(
    monos,
    sp::Union{SignSymmetry,Monomial},
    gram_basis::MB.SubBasis{MB.Monomial},
)
    return MB.SubBasis{MB.Monomial}.(sparsity(monos, sp, gram_basis.monomials))
end
function sparsity(
    poly::MP.AbstractPolynomial,
    sp::Union{SignSymmetry,Monomial},
    certificate::SumOfSquares.Certificate.AbstractIdealCertificate,
)
    return sparsity(
        MP.monomials(poly),
        sp,
        SumOfSquares.Certificate.gram_basis(certificate, poly),
    )
end
function sparsity(v::SumOfSquares.Certificate.WithVariables, sp, certificate)
    return sparsity(v.inner, sp, certificate)
end
function SumOfSquares.Certificate.gram_basis(certificate::Ideal, poly)
    return sparsity(poly, certificate.sparsity, certificate.certificate)
end
function SumOfSquares.Certificate.gram_basis_type(
    ::Type{Ideal{S,C}},
    ::Type{M},
) where {S,C,M}
    return SumOfSquares.Certificate.gram_basis_type(C, M)
end
function SumOfSquares.Certificate.reduced_polynomial(
    certificate::Ideal,
    poly,
    domain,
)
    return SumOfSquares.Certificate.reduced_polynomial(
        certificate.certificate,
        poly,
        domain,
    )
end
function SumOfSquares.Certificate.reduced_basis(certificate::Ideal, basis, domain)
    return SumOfSquares.Certificate.reduced_basis(certificate.certificate, basis, domain)
end
function MA.promote_operation(::typeof(SumOfSquares.Certificate.reduced_basis), ::Type{Ideal{S,C}}, ::Type{B}, ::Type{D}) where {S,C,B,D}
    return MA.promote_operation(SumOfSquares.Certificate.reduced_basis, C, B, D)
end
function SumOfSquares.Certificate.cone(certificate::Ideal)
    return SumOfSquares.Certificate.cone(certificate.certificate)
end
function SumOfSquares.matrix_cone_type(::Type{Ideal{S,C}}) where {S,C}
    return SumOfSquares.matrix_cone_type(C)
end

function SumOfSquares.Certificate.zero_basis(certificate::Ideal)
    return SumOfSquares.Certificate.zero_basis(certificate.certificate)
end
function SumOfSquares.Certificate.zero_basis_type(
    ::Type{Ideal{S,C}},
) where {S,C}
    return SumOfSquares.Certificate.zero_basis_type(C)
end
