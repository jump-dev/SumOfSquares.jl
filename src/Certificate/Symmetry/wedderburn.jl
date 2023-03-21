function SymbolicWedderburn.decompose(
    k::MP.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if features(hom) == monomials

    indcs = [hom[mono] for mono in MP.monomials(k)]
    coeffs = MP.coefficients(k)

    return indcs, coeffs
end

function SymbolicWedderburn.ExtensionHomomorphism(
    action::SymbolicWedderburn.Action,
    basis::MB.MonomialBasis,
)
    monos = collect(basis.monomials)
    return SymbolicWedderburn.ExtensionHomomorphism(Int, action, monos)
end

struct VariablePermutation <: SymbolicWedderburn.ByPermutations end
_map_idx(f, v::AbstractVector) = map(f, eachindex(v))
_tuple_map_idx(f, ::Tuple{}, i) = tuple()
function _tuple_map_idx(f, v::Tuple, i)
    return (f(i), _tuple_map_idx(f, Base.tail(v), i + 1)...)
end
_map_idx(f, v::Tuple) = _tuple_map_idx(f, v, 1)
function SymbolicWedderburn.action(
    ::VariablePermutation,
    p,
    mono::MP.AbstractMonomial,
)
    v = MP.variables(mono)
    return MP.substitute(MP.Eval(), mono, v => _map_idx(i -> v[i^p], v))
end
abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el,
    term::MP.AbstractTerm,
)
    return MP.coefficient(term) *
           SymbolicWedderburn.action(a, el, MP.monomial(term))
end
function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el,
    poly::MP.AbstractPolynomial,
)
    return sum([
        SymbolicWedderburn.action(a, el, term) for term in MP.terms(poly)
    ])
end

# TODO Move it to MultivariateBases
function MP.polynomialtype(
    ::Type{<:MB.AbstractPolynomialVectorBasis{PT}},
    T::Type,
) where {PT}
    C = MP.coefficienttype(PT)
    U = MA.promote_operation(*, C, T)
    V = MA.promote_operation(+, U, U)
    return MP.polynomialtype(PT, V)
end

"""
    struct Symmetry.Ideal{C,GT,AT<:SymbolicWedderburn.Action} <: SumOfSquares.Certificate.AbstractIdealCertificate
        pattern::Symmetry.Pattern{GT,AT}
        certificate::C
    end

Same certificate as `certificate` except that the Sum-of-Squares polynomial `σ`
is modelled as a sum of Sum-of-Squares polynomials with smaller bases
using the Symbolic Wedderburn decomposition of the symmetry pattern specified
by `pattern`.
"""
struct Ideal{C,GT,AT<:SymbolicWedderburn.Action} <:
       SumOfSquares.Certificate.AbstractIdealCertificate
    pattern::Pattern{GT,AT}
    certificate::C
end
function SumOfSquares.Certificate.cone(certificate::Ideal)
    return SumOfSquares.Certificate.cone(certificate.certificate)
end
function SumOfSquares.matrix_cone_type(::Type{<:Ideal{C}}) where {C}
    return SumOfSquares.matrix_cone_type(C)
end
function SumOfSquares.Certificate.gram_basis_type(::Type{<:Ideal})
    return Vector{Vector{MB.FixedPolynomialBasis}}
end
SumOfSquares.Certificate.zero_basis_type(::Type{<:Ideal}) = MB.MonomialBasis
SumOfSquares.Certificate.zero_basis(::Ideal) = MB.MonomialBasis
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

function matrix_reps(pattern, R, basis, ::Type{T}, form) where {T}
    polys = R * basis.monomials
    return map(SymbolicWedderburn.gens(pattern.group)) do g
        S = Matrix{T}(undef, length(polys), length(polys))
        for i in eachindex(polys)
            p = polys[i]
            q = SymbolicWedderburn.action(pattern.action, g, p)
            coefs = MP.coefficients(q, basis.monomials)
            S[:, i] = _linsolve(R, coefs, form)
        end
        return S
    end
end

function SumOfSquares.Certificate.gram_basis(cert::Ideal, poly)
    basis = SumOfSquares.Certificate.gram_basis(cert.certificate, poly)
    T = SumOfSquares._complex(
        Float64,
        SumOfSquares.matrix_cone_type(typeof(cert)),
    )
    return _gram_basis(cert.pattern, basis, T)
end

function _gram_basis(pattern::Pattern, basis, ::Type{T}) where {T}
    # We set `semisimple=true` as we don't support simple yet since it would not give all the simple components but only one of them.
    summands = SymbolicWedderburn.symmetry_adapted_basis(
        T,
        pattern.group,
        pattern.action,
        basis,
        semisimple = true,
    )
    # We have a new basis `b = vcat(R * basis.monomials for R in summands)``.
    # SymbolicWedderburn guarantees that the invariant subspace spanned by the
    # polynomials of the vector `R * basis.monomials` is invariant under the
    # action of the group. That is, the matrix representation `ρ(g)` induced by
    # this basis is block-diagonal with one block for each summand.
    # That block is the matrix `S` computed below.
    # So an invariant solution `b'*Q*b` satisfies `Diagonal(S' for S in ...) * Q * Diagonal(S for S in ...) = Q`.
    # Or in equivalently: `Q * Diagonal(S for S in ...) = Diagonal(inv(S') for S in ...) * Q`.
    form = if T <: Union{AbstractFloat,Complex{<:AbstractFloat}}
        _OrthogonalMatrix()
    else
        _RowEchelonMatrix()
    end
    return map(summands) do summand
        R = SymbolicWedderburn.image_basis(summand)
        m = SymbolicWedderburn.multiplicity(summand)
        N = size(R, 1)
        d = SymbolicWedderburn.degree(summand)
        S = matrix_reps(pattern, R, basis, T, form)
        #S = matrix_reps(pattern, R, basis, T, _RowEchelonMatrix())
        decomose_semisimple = d > 1
        if decomose_semisimple
            # If it's not orthogonal, how can we conclude that we can still use the semisimple summands block-decomposition ?
            # In Example 1.7.2 of Sagan's book, he uses Corollary 1.6.6 which requires that `X` and `Y` are irreducible
            # (where `X` and `Y` are here the `S` corresponding to two different summands).
            # Here, given semisimple representations `X` and `Y`, they are not irreducible if `m > 1`.
            # Furthermore, as they are not orthogonal, we have something like `T * X = inv(Y') * T`, so how can we know that
            # `X` and `inv(Y')` are not equivalent (to exclude the case 1. of Corollary 1.6.6) ?
            if !all(is_orthogonal, S)
                R = orthogonalize(R)
                S = matrix_reps(pattern, R, basis, T, _OrthogonalMatrix())
                for i in axes(R, 1)
                    R[i, :] = LinearAlgebra.normalize(R[i, :])
                end
                S = matrix_reps(pattern, R, basis, T, _OrthogonalMatrix())
                if !all(is_orthogonal, S)
                    error(
                        "The matrix representation induced from the action on the polynomial basis is not orthogonal.",
                    )
                    # We would like to just throw this warning and just not decompose the semisimple summand but
                    # as explained in the comment above, it's not even clear that the diagonalization induced by the simple summands is correct.
                    #@warn("The matrix representation induced from the action on the polynomial basis is not orthogonal. The $(m * d)-dimensional semisimple summand can be decomposed onto $m simple summands of degree $d so that the $(m * d) x $(m * d) diagonal block is reduced to $d identical copied of a single $m x $m diagonal block. However, as the action is not orthogonal, this decomposition will not happen.")
                    #decomose_semisimple = false
                end
            end
        end
        F = convert(Matrix{T}, R)
        if d > 1
            if m > 1
                U = ordered_block_diag(S, d)
            else
                U = Matrix{T}(LinearAlgebra.I, N, N)
            end
            if U === nothing
                error(
                    "Could not simultaneously block-diagonalize into $m identical $dx$d blocks",
                )
            end
            # From Example 1.7.3 of
            # Sagan, The symmetric group, Springer Science & Business Media, 2001
            # we know that there exists `C` such that `Q = kron(C, I)` if we use
            # `(U[1:d] * F)' * basis.monomials`, `(U[d+1:2d] * F)' * basis.monomials`, ...
            # where `C` are some complex numbers as they are eigenvalues (see Corollary 1.6.8).
            # As `Q` is symmetric, we know the eigenvalues are real so we can take `C` real as well.
            # Moreover, `Q = kron(C, I)` is not block diagonal but we can get a block-diagonal
            # `Q = kron(I, Q)` by permuting the rows and columns:
            # `(U[1:d:(1+d*(m-1))] * F)' * basis.monomials`, `(U[2:d:(2+d*(m-1))] * F)' * basis.monomials`, ...
            map(1:d) do i
                return MB.FixedPolynomialBasis(
                    (transpose(U[:, i:d:(i+d*(m-1))]) * F) * basis.monomials,
                )
            end
        else
            F = convert(Matrix{T}, R)
            [MB.FixedPolynomialBasis(F * basis.monomials)]
        end
    end
end
