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
    mono_to_index = Dict(monos[i] => i for i in eachindex(monos))
    return SymbolicWedderburn.ExtensionHomomorphism(action, monos, mono_to_index)
end

struct VariablePermutation <: SymbolicWedderburn.ByPermutations end
_map_idx(f, v::AbstractVector) = map(f, eachindex(v))
_tuple_map_idx(f, ::Tuple{}, i) = tuple()
_tuple_map_idx(f, v::Tuple, i) = (f(i), _tuple_map_idx(f, Base.tail(v), i + 1)...)
_map_idx(f, v::Tuple) = _tuple_map_idx(f, v, 1)
function SymbolicWedderburn.action(::VariablePermutation, p, mono::MP.AbstractMonomial)
    v = MP.variables(mono)
    MP.substitute(MP.Eval(), mono, v => _map_idx(i -> v[i^p], v))
end
abstract type OnMonomials <: SymbolicWedderburn.ByLinearTransformation end

function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el,
    term::MP.AbstractTerm,
)
    return MP.coefficient(term) * SymbolicWedderburn.action(a, el, MP.monomial(term))
end
function SymbolicWedderburn.action(
    a::Union{VariablePermutation,OnMonomials},
    el,
    poly::MP.AbstractPolynomial,
)
    return sum([SymbolicWedderburn.action(a, el, term) for term in MP.terms(poly)])
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

struct Ideal{C,GT,AT<:SymbolicWedderburn.Action} <: SumOfSquares.Certificate.AbstractIdealCertificate
    pattern::Pattern{GT,AT}
    certificate::C
end
SumOfSquares.Certificate.get(certificate::Ideal, attr::SumOfSquares.Certificate.Cone) = SumOfSquares.Certificate.get(certificate.certificate, attr)
function SumOfSquares.matrix_cone_type(::Type{<:Ideal{C}}) where {C}
    return SumOfSquares.matrix_cone_type(C)
end
SumOfSquares.Certificate.get(::Type{<:Ideal}, ::SumOfSquares.Certificate.GramBasisType) = Vector{Vector{MB.FixedPolynomialBasis}}
SumOfSquares.Certificate.zero_basis_type(::Type{<:Ideal}) = MB.MonomialBasis
SumOfSquares.Certificate.zero_basis(::Ideal) = MB.MonomialBasis
function SumOfSquares.Certificate.get(certificate::Ideal, attr::SumOfSquares.Certificate.ReducedPolynomial, poly, domain)
    return SumOfSquares.Certificate.get(certificate.certificate, attr, poly, domain)
end

function matrix_reps(group_els, action, R, basis, ::Type{T}) where {T}
    polys = R * basis
    return map(group_els) do g
        S = Matrix{T}(undef, length(polys), length(polys))
        for i in eachindex(polys)
            p = polys[i]
            q = SymbolicWedderburn.action(action, g, p)
            coefs = MP.coefficients(q, basis)
            col = row_echelon_linsolve(R, coefs)
            S[:, i] = col
        end
        return S
    end
end

function check_orthogonal(S, tol=1e-6)
    n = LinearAlgebra.checksquare(S)
    for i in 1:n
        for j in 1:n
            if i != j
                s = LinearAlgebra.dot(S[:, i], S[:, j])
                if abs(s) > tol
                    @warn("Matrix not orthogonal, scalar product between column $i and $j is $s. Expect symmetry reduction to be conservative.")
                end
            end
        end
    end
end

function SumOfSquares.Certificate.get(cert::Ideal, attr::SumOfSquares.Certificate.GramBasis, poly)
    basis = SumOfSquares.Certificate.get(cert.certificate, attr, poly)
    T = SumOfSquares._complex(Float64, SumOfSquares.matrix_cone_type(typeof(cert)))
    summands = SymbolicWedderburn.symmetry_adapted_basis(T, cert.pattern.group, basis, cert.pattern.action)
    return map(summands) do summand
        R = SymbolicWedderburn.basis(summand)
        m = SymbolicWedderburn.multiplicity(summand)
        F = convert(Matrix{T}, R)
        N = size(R, 1)
        d = SymbolicWedderburn.degree(summand)
        S = matrix_reps(SymbolicWedderburn.gens(cert.pattern.group), cert.pattern.action, R, basis.monomials, T)
        check_orthogonal.(S)
        if d > 1
            if m > 1
                U = ordered_block_diag(S, d)
            else
                U = Matrix{T}(LinearAlgebra.I, N, N)
            end
            if U === nothing
                error("Could not simultaneously block-diagonalize into $m identical $(d)x$d blocks")
            end
            # From Example 1.7.3 of
            # Sagan, The symmetric group, Springer Science & Business Media, 2001
            # we know that there exists `C` such that `Q = kron(C, I)` if we use
            # `(U[1:d] * F)' * basis.monomials`, `(U[d+1:2d] * F)' * basis.monomials`, ...
            # where `C` are some complex numbers as they are eigenvalues (see Corollary 1.6.8).
            # As `Q` is symmetric, we now the eigenvalues are real so we can take `C` real as well.
            # Moreover, `Q = kron(C, I)` is not block diagonal but we can get a block-diagonal
            # `Q = kron(I, Q)` by permuting the rows and columns:
            # `(U[1:d:(1+d*(m-1))] * F)' * basis.monomials`, `(U[2:d:(2+d*(m-1))] * F)' * basis.monomials`, ...
            map(1:d) do i
                MB.FixedPolynomialBasis(
                    (transpose(U[:, i:d:(i + d * (m - 1))]) * F) * basis.monomials,
                )
            end
        else
            [MB.FixedPolynomialBasis(F * basis.monomials)]
        end
    end
end
