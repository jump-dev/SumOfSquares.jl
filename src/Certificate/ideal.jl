######################
# Ideal certificates #
######################

abstract type AbstractIdealCertificate <: AbstractCertificate end

struct _NonZero <: Number end
Base.iszero(::_NonZero) = false
Base.convert(::Type{_NonZero}, ::Number) = _NonZero()
Base.:*(a::_NonZero, ::Number) = a
Base.:*(::Number, a::_NonZero) = a
Base.:*(::_NonZero, a::_NonZero) = a
Base.:+(a::_NonZero, ::Number) = a
Base.:+(::Number, a::_NonZero) = a
Base.:+(::_NonZero, a::_NonZero) = a

function _combine_with_gram(
    basis::MB.SubBasis{B,M},
    gram_bases::AbstractVector{<:MB.SubBasis},
    weights,
) where {B,M}
    full = MB.FullBasis{B,M}()
    mstr = SA.mstructure(full)
    p = MP.polynomial(fill(_NonZero(), length(basis)), basis.monomials)
    for (gram, weight) in zip(gram_bases, weights)
        w = MB.convert_basis(full, weight)
        for j in eachindex(gram)
            for i in 1:j
                s = MB.convert_basis(full, gram[i] * gram[j])
                for w_mono in SA.supp(w)
                    for s_mono in SA.supp(s)
                        MA.operate!(
                            SA.UnsafeAddMul(mstr),
                            p,
                            w_mono.monomial,
                            s_mono.monomial,
                        )
                    end
                end
            end
        end
    end
    MA.operate!(SA.canonical, p)
    return MB.SubBasis{MB.Monomial}(MP.monomials(p))
end

_reduce_with_domain(basis::MB.SubBasis{MB.Monomial}, ::FullSpace) = basis

function _reduce_with_domain(basis::MB.SubBasis{MB.Monomial}, domain)
    I = ideal(domain)
    # set of standard monomials that are hit
    standard = Set{eltype(basis.monomials)}()
    for mono in basis.monomials
        r = rem(mono, I)
        union!(standard, MP.monomials(r))
    end
    return MB.QuotientBasis(
        MB.SubBasis{MB.Monomial}(MP.monomial_vector(collect(standard))),
        I,
    )
end

function reduced_basis(
    ::AbstractIdealCertificate,
    basis,
    domain,
    gram_bases,
    weights,
)
    return _reduce_with_domain(
        _combine_with_gram(basis, gram_bases, weights),
        domain,
    )
end

abstract type SimpleIdealCertificate{C,B} <: AbstractIdealCertificate end

reduced_polynomial(::SimpleIdealCertificate, poly, domain) = poly

cone(certificate::SimpleIdealCertificate) = certificate.cone
function SumOfSquares.matrix_cone_type(
    ::Type{<:SimpleIdealCertificate{CT}},
) where {CT}
    return SumOfSquares.matrix_cone_type(CT)
end

# TODO return something else when `PolyJuMP` support other bases.
zero_basis(certificate::SimpleIdealCertificate) = certificate.basis # FIXME not used yet
function zero_basis_type(::Type{<:SimpleIdealCertificate{C,B}}) where {C,B}
    return B
end

"""
    struct MaxDegree{C<:SumOfSquares.SOSLikeCone,B<:SA.AbstractBasis} <: SimpleIdealCertificate{C,B}
        cone::C
        basis::B
        maxdegree::Int
    end

The `MaxDegree` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with a basis of type `basis` such that
the degree of `σ(x)` does not exceed `maxdegree`.
"""
struct MaxDegree{C<:SumOfSquares.SOSLikeCone,B<:SA.AbstractBasis} <:
       SimpleIdealCertificate{C,B}
    cone::C
    basis::B
    maxdegree::Int
end
function gram_basis(certificate::MaxDegree, poly)
    return maxdegree_gram_basis(
        certificate.basis,
        MP.variables(poly),
        certificate.maxdegree,
    )
end
function gram_basis_type(::Type{MaxDegree{C,B}}, ::Type{M}) where {C,B,M} # TODO remove `M` it's not needed anymore
    return MB.explicit_basis_type(B)
end

"""
    struct FixedBasis{C<:SumOfSquares.SOSLikeCone,B<:SA.ExplicitBasis} <: SimpleIdealCertificate{C,B}
        cone::CT
        basis::BT
    end

The `FixedBasis` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with basis `basis`.
"""
struct FixedBasis{C<:SumOfSquares.SOSLikeCone,B<:SA.ExplicitBasis} <:
       SimpleIdealCertificate{C,B}
    cone::C
    basis::B
end
function gram_basis(certificate::FixedBasis, _)
    return certificate.basis
end
gram_basis_type(::Type{FixedBasis{C,B}}, ::Type) where {C,B} = B

"""
    struct Newton{
        C<:SumOfSquares.SOSLikeCone,
        B<:SA.AbstractBasis,
        N<:AbstractNewtonPolytopeApproximation,
    } <: SimpleIdealCertificate{C,B}
        cone::C
        basis::B
        newton::N
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
struct Newton{
    C<:SumOfSquares.SOSLikeCone,
    B<:SA.AbstractBasis,
    N<:AbstractNewtonPolytopeApproximation,
} <: SimpleIdealCertificate{C,B}
    cone::C
    basis::B
    newton::N
end

function Newton(cone, basis, variable_groups::Tuple)
    return Newton(
        cone,
        basis,
        NewtonFilter(NewtonDegreeBounds(variable_groups)),
    )
end

function gram_basis(certificate::Newton, basis)
    return MB.basis_covering_monomials(
        certificate.basis,
        monomials_half_newton_polytope(SA.basis(basis), certificate.newton),
    )
end

function gram_basis_type(::Type{<:Newton{C,B}}, ::Type{M}) where {C,B,M}
    return MB.explicit_basis_type(B)
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

function _rem(coeffs, basis::MB.SubBasis{MB.Monomial}, I)
    poly = MP.polynomial(coeffs, basis.monomials)
    return convert(typeof(poly), rem(poly, I))
end

function reduced_polynomial(
    ::Remainder,
    a::SA.AlgebraElement,
    domain,
)
    r = _rem(SA.coeffs(a), SA.basis(a), ideal(domain))
    return MB.algebra_element(
        MP.coefficients(r),
        MB.SubBasis{MB.Monomial}(MP.monomials(r)),
    )
end

function gram_basis(certificate::Remainder, poly)
    return gram_basis(certificate.gram_certificate, poly)
end

function gram_basis_type(::Type{Remainder{GCT}}, ::Type{M}) where {GCT,M}
    return gram_basis_type(GCT, M)
end

cone(certificate::Remainder) = cone(certificate.gram_certificate)
function SumOfSquares.matrix_cone_type(::Type{Remainder{GCT}}) where {GCT}
    return SumOfSquares.matrix_cone_type(GCT)
end
zero_basis(certificate::Remainder) = zero_basis(certificate.gram_certificate)
zero_basis_type(::Type{Remainder{GCT}}) where {GCT} = zero_basis_type(GCT)

function _quotient_basis_type(
    ::Type{B},
    ::Type{D},
) where {T,I,B<:SA.AbstractBasis{T,I},D}
    return MB.QuotientBasis{
        T,
        I,
        B,
        MA.promote_operation(SemialgebraicSets.ideal, D),
    }
end

function MA.promote_operation(
    ::typeof(reduced_basis),
    ::Type{<:Union{SimpleIdealCertificate,Remainder}},
    ::Type{B},
    ::Type{SemialgebraicSets.FullSpace},
    ::Type,
    ::Type,
) where {B}
    return B
end

function MA.promote_operation(
    ::typeof(reduced_basis),
    ::Type{<:Union{SimpleIdealCertificate,Remainder}},
    ::Type{B},
    ::Type{D},
    ::Type,
    ::Type,
) where {B,D}
    return _quotient_basis_type(B, D)
end
