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
    basis::MB.SubBasis{B},
    gram_bases::AbstractVector{<:SA.ExplicitBasis},
    weights,
) where {B}
    for gram_basis in gram_bases
        if parent(basis) != parent(gram_basis)
            error(
                "Bases $(parent(basis)) and $(parent(gram_basis)) are incompatible",
            )
        end
    end
    p = zero(_NonZero, MB.algebra(parent(basis)))
    cache = zero(_NonZero, MB.algebra(parent(basis)))
    for mono in basis
        MA.operate!(
            SA.UnsafeAddMul(*),
            p,
            _term_constant_monomial(_NonZero(), mono),
            MB.algebra_element(mono),
        )
    end
    for (gram, weight) in zip(gram_bases, weights)
        if isone(weight)
            MA.operate!(
                SA.UnsafeAdd(),
                p,
                SA.QuadraticForm(
                    GramMatrix{_NonZero}((_, _) -> _NonZero(), gram),
                ),
            )
        else
            MA.operate_to!(
                cache,
                +,
                GramMatrix{_NonZero}((_, _) -> _NonZero(), gram),
            )
            MA.operate!(SA.UnsafeAddMul(*), p, cache, weight)
        end
    end
    MA.operate!(SA.canonical, SA.coeffs(p))
    return SA.sub_basis(parent(basis), keys(SA.coeffs(p)))
end

function _reduce_with_domain(basis::MB.SubBasis, zero_basis, ::FullSpace)
    return MB.explicit_basis_covering(zero_basis, basis)
end

function _reduce_with_domain(basis, zero_basis, domain)
    return __reduce_with_domain(basis, zero_basis, domain)
end

function __reduce_with_domain(_, _, _)
    return error("Only Monomial basis support with an equalities in domain")
end

function __reduce_with_domain(
    basis::MB.SubBasis{MB.Monomial},
    zero_basis::MB.FullBasis{MB.Monomial},
    domain,
)
    return __reduce_with_domain(basis, zero_basis, ideal(domain))
end

function __reduce_with_domain(
    basis::MB.SubBasis{MB.Monomial},
    ::MB.FullBasis{MB.Monomial},
    I::SemialgebraicSets.AbstractPolynomialIdeal,
)
    T = float(MP.coefficient_type(typeof(I)))
    return MB.QuotientBasis(basis, _remainder_matrix(basis, I, T))
end

# The Gröbner basis of the ideal of an `AlgebraicSet` may not be computed and
# computing it may be expensive. In that case, we consider instead the
# condition that the polynomial is in the row space of the Macaulay matrix
# of the equalities of degree `maxdegree(basis)`. This corresponds to
# searching for multipliers `λ_i` of the equalities `h_i` such that the degree
# of `λ_i * h_i` is at most `maxdegree(basis)`.
# If the Gröbner basis is already computed, we use the remainder modulo the
# ideal instead. In both cases, the result is a `QuotientBasis` with a sparse
# matrix so that the type of the zero basis can be inferred from the type of
# the domain.
function __reduce_with_domain(
    basis::MB.SubBasis{MB.Monomial},
    ::MB.FullBasis{MB.Monomial},
    domain::AlgebraicSet{T},
) where {T}
    I = ideal(domain)
    if SemialgebraicSets.is_gröbner_basis(I)
        W = _remainder_matrix(basis, I, float(T))
    else
        W = _macaulay_zero_matrix(basis, equalities(domain), float(T))
    end
    return MB.QuotientBasis(basis, W)
end

# Sparse matrix `R` such that `R[i, j]` is the coefficient of the `i`th
# standard monomial in the remainder of the `j`th monomial of `basis` modulo
# `I`. Its number of nonzeros is the total number of terms of the remainders.
function _remainder_matrix(basis, I, ::Type{T}) where {T}
    if isempty(basis)
        return SparseArrays.spzeros(T, 0, 0)
    end
    rems = [rem(mono, I) for mono in MB.keys_as_monomials(basis)]
    standard = MP.monomial_vector(reduce(vcat, MP.monomials.(rems)))
    rows = Int[]
    cols = Int[]
    vals = T[]
    for (j, r) in enumerate(rems)
        for t in MP.terms(r)
            push!(rows, searchsortedfirst(standard, MP.monomial(t)))
            push!(cols, j)
            push!(vals, MP.coefficient(t))
        end
    end
    return SparseArrays.sparse(
        rows,
        cols,
        vals,
        length(standard),
        length(basis),
    )
end

# Sparse matrix `W` such that a polynomial `t` with support in `basis` is in
# the row space of the Macaulay matrix of degree `d` of `equalities` if and
# only if `W * coefficients(t, basis) = 0` where `d = maxdegree(basis)`.
function _macaulay_zero_matrix(basis, equalities, ::Type{T}) where {T}
    n = length(basis)
    if isempty(basis)
        return SparseArrays.spzeros(T, 0, 0)
    end
    monos = MB.keys_as_monomials(basis)
    # The shifts should be in all the variables and not only in the variables
    # of the equalities
    vars =
        MP.variables(MP.monomial_vector([monos; MP.monomials.(equalities)...]))
    M = MB.macaulay(equalities, MP.maxdegree(basis), vars)
    null = MM.MacaulayNullspace(M)
    columns = MP.variables(null.basis)
    rows = Int[]
    found = Int[]
    others = Int[]
    for (i, mono) in enumerate(monos)
        j = SA.key_index(null.basis, MP.exponents(mono, columns))
        if isnothing(j)
            # The monomial is not in any row of the Macaulay matrix so it is
            # independent of the other ones
            push!(others, i)
        else
            push!(rows, j)
            push!(found, i)
        end
    end
    Z = null.matrix[rows, :]
    if isempty(Z)
        r = 0
    else
        F = LinearAlgebra.svd(Z)
        rank_check = MM.LeadingRelativeRankTol(Base.rtoldefault(T))
        r = MM.rank_from_singular_values(F.S, rank_check)
    end
    if r == length(found)
        # No monomial of `basis` is dependent so the polynomial must be zero
        return SparseArrays.sparse(one(T) * LinearAlgebra.I, n, n)
    end
    # The equations `transpose(F.U[:, 1:r]) * coefficients(t, basis) = 0` are
    # orthonormal but they are dense which is an issue for sparse solvers.
    # As any invertible recombination of the equations is equivalent, we
    # take the reduced row echelon form instead. The pivots are the smallest
    # independent monomials so if they coincide with the standard monomials
    # of a Gröbner basis, it coincides with the matrix of `_remainder_matrix`.
    E = Matrix(transpose(F.U[:, 1:r]))
    tol = Base.rtoldefault(T)
    RowEchelon.rref!(E, tol)
    # Remove the numerical noise so that the zeros are skipped
    E[abs.(E) .< tol] .= zero(T)
    # The pivot of each row of `W`, we sort the rows according to their
    # pivots so that the equations are in the order of `basis`
    pivots = [
        [found[findfirst(!iszero, view(E, k, :))] for k in 1:r]
        others
    ]
    W = zeros(T, length(pivots), n)
    W[1:r, found] = E
    for (k, i) in enumerate(others)
        W[r+k, i] = one(T)
    end
    return SparseArrays.sparse(W[sortperm(pivots), :])
end

function zero_basis(
    cert::AbstractIdealCertificate,
    basis,
    domain,
    gram_bases,
    weights,
)
    return _reduce_with_domain(
        _combine_with_gram(basis, gram_bases, weights),
        _zero_basis(cert),
        domain,
    )
end

abstract type SimpleIdealCertificate{C,G,Z} <: AbstractIdealCertificate end

reduced_polynomial(::SimpleIdealCertificate, poly, domain) = poly

cone(certificate::SimpleIdealCertificate) = certificate.cone
function SumOfSquares.matrix_cone_type(
    ::Type{<:SimpleIdealCertificate{CT}},
) where {CT}
    return SumOfSquares.matrix_cone_type(CT)
end

function MA.promote_operation(
    ::typeof(gram_basis),
    ::Type{<:SimpleIdealCertificate{C,G}},
) where {C,G}
    return MB.explicit_basis_type(G)
end

"""
    struct MaxDegree{C<:SumOfSquares.SOSLikeCone,G<:SA.AbstractBasis,Z<:SA.AbstractBasis} <:
        SimpleIdealCertificate{C,G,Z}
        cone::C
        gram_basis::G
        zero_basis::Z
        maxdegree::Int
    end

The `MaxDegree` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with a basis of type `basis` such that
the degree of `σ(x)` does not exceed `maxdegree`.
"""
struct MaxDegree{
    C<:SumOfSquares.SOSLikeCone,
    G<:SA.AbstractBasis,
    Z<:SA.AbstractBasis,
} <: SimpleIdealCertificate{C,G,Z}
    cone::C
    gram_basis::G
    zero_basis::Z
    maxdegree::Int
end

MB.implicit_basis(certificate::MaxDegree) = certificate.gram_basis

function gram_basis(certificate::MaxDegree, poly)
    return maxdegree_gram_basis(
        certificate.gram_basis,
        MP.variables(poly),
        certificate.maxdegree,
    )
end

"""
    struct FixedBasis{C<:SumOfSquares.SOSLikeCone,G<:SA.ExplicitBasis,Z<:SA.AbstractBasis} <:
        SimpleIdealCertificate{C,G,Z}
        cone::C
        gram_basis::G
        zero_basis::Z
    end

The `FixedBasis` certificate ensures the nonnegativity of `p(x)` for all `x` such that
`h_i(x) = 0` by exhibiting a Sum-of-Squares polynomials `σ(x)`
such that `p(x) - σ(x)` is guaranteed to be zero for all `x`
such that `h_i(x) = 0`.
The polynomial `σ(x)` is search over `cone` with basis `basis`.
"""
struct FixedBasis{
    C<:SumOfSquares.SOSLikeCone,
    G<:SA.ExplicitBasis,
    Z<:SA.AbstractBasis,
} <: SimpleIdealCertificate{C,G,Z}
    cone::C
    gram_basis::G
    zero_basis::Z
end

function MB.implicit_basis(certificate::FixedBasis)
    return MB.implicit_basis(certificate.gram_basis)
end

function gram_basis(certificate::FixedBasis, _)
    return certificate.gram_basis
end

function MA.promote_operation(
    ::typeof(gram_basis),
    ::Type{<:FixedBasis{C,G}},
) where {C,G}
    return G
end

"""
    struct Newton{
        C<:SumOfSquares.SOSLikeCone,
        G<:SA.AbstractBasis,
        Z<:SA.AbstractBasis,
        N<:AbstractNewtonPolytopeApproximation,
    } <: SimpleIdealCertificate{C,G,Z}
        cone::C
        gram_basis::G
        zero_basis::Z
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
    G<:SA.AbstractBasis,
    Z<:SA.AbstractBasis,
    N<:AbstractNewtonPolytopeApproximation,
} <: SimpleIdealCertificate{C,G,Z}
    cone::C
    gram_basis::G
    zero_basis::Z
    newton::N
end

function Newton(cone, gram_basis, zero_basis, variable_groups::Tuple)
    return Newton(
        cone,
        gram_basis,
        zero_basis,
        NewtonFilter(NewtonDegreeBounds(variable_groups)),
    )
end

MB.implicit_basis(certificate::Newton) = certificate.gram_basis

function gram_basis(certificate::Newton, poly)
    basis = half_newton_polytope(
        _algebra_element(poly),
        MP.variables(poly),
        certificate.newton,
    )
    return MB.explicit_basis_covering(certificate.gram_basis, basis)
end

"""
    struct Remainder{C<:AbstractIdealCertificate} <: AbstractIdealCertificate
        gram_certificate::C
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
struct Remainder{C<:AbstractIdealCertificate} <: AbstractIdealCertificate
    gram_certificate::C
end

function _rem(coeffs, basis::MB.FullBasis{MB.Monomial}, I)
    poly = MP.polynomial(MB.algebra_element(coeffs, basis))
    r = convert(typeof(poly), rem(poly, I))
    # `FixedVariablesSet`, use substitutions which makes us
    # loose variables
    s, _ = SA.promote_bases(r, poly)
    return MB.algebra_element(SA.coeffs(s), basis)
end

function reduced_polynomial(::Remainder, a::SA.AlgebraElement, domain)
    return _rem(SA.coeffs(a), SA.basis(a), ideal(domain))
end

function MB.implicit_basis(certificate::Remainder)
    return MB.implicit_basis(certificate.gram_certificate)
end

function gram_basis(certificate::Remainder, poly)
    return gram_basis(certificate.gram_certificate, poly)
end

function MA.promote_operation(
    ::typeof(gram_basis),
    ::Type{Remainder{C}},
) where {C}
    return MA.promote_operation(gram_basis, C)
end

cone(certificate::Remainder) = cone(certificate.gram_certificate)

function SumOfSquares.matrix_cone_type(::Type{Remainder{GCT}}) where {GCT}
    return SumOfSquares.matrix_cone_type(GCT)
end

_zero_basis(c::SimpleIdealCertificate) = c.zero_basis

function _zero_basis_type(::Type{<:SimpleIdealCertificate{C,G,Z}}) where {C,G,Z}
    return Z
end

_zero_basis(c::Remainder) = _zero_basis(c.gram_certificate)

# The polynomial was already reduced with `rem` by `reduced_polynomial` so the
# Gröbner basis is computed and we use the remainder matrix.
function zero_basis(cert::Remainder, basis, domain, gram_bases, weights)
    return _reduce_with_domain(
        _combine_with_gram(basis, gram_bases, weights),
        _zero_basis(cert),
        _remainder_domain(domain),
    )
end

_remainder_domain(domain) = ideal(domain)
_remainder_domain(domain::FullSpace) = domain

function _zero_basis_type(::Type{Remainder{C}}) where {C}
    return _zero_basis_type(C)
end

function MA.promote_operation(
    ::typeof(zero_basis),
    C::Type{<:Union{SimpleIdealCertificate,Remainder}},
    ::Type,
    ::Type{SemialgebraicSets.FullSpace},
    ::Type,
    ::Type,
)
    return MB.explicit_basis_type(_zero_basis_type(C))
end

function MA.promote_operation(
    ::typeof(zero_basis),
    ::Type{<:Union{SimpleIdealCertificate,Remainder}},
    ::Type{B},
    ::Type{D},
    ::Type,
    ::Type,
) where {B,D}
    return _quotient_basis_type(B, float(MP.coefficient_type(D)))
end

function _quotient_basis_type(
    ::Type{B},
    ::Type{T},
) where {E,I,B<:SA.AbstractBasis{E,I},T}
    return MB.QuotientBasis{E,I,B,SparseArrays.SparseMatrixCSC{T,Int}}
end
