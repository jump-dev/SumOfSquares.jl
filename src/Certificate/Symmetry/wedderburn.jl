function SymbolicWedderburn.decompose(
    p::MB.Polynomial,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    return [hom[p]], [1]
end

function SymbolicWedderburn.decompose(
    p::SA.AlgebraElement,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    return [hom[k] for k in SA.supp(p)],
    [v for (_, v) in SA.nonzero_pairs(SA.coeffs(p))]
end

function SymbolicWedderburn.decompose(
    k::MP.AbstractPolynomialLike,
    hom::SymbolicWedderburn.InducedActionHomomorphism,
)
    # correct only if features(hom) == monomials

    indcs = [hom[mono] for mono in MP.monomials(k)]
    coeffs = MP.coefficients(k)

    return indcs, coeffs
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
    p::MB.Polynomial{MB.Monomial},
)
    res = SymbolicWedderburn.action(a, el, MP.monomial(p))
    if res isa MP.AbstractMonomial
        return MB.Polynomial{MB.Monomial}(res)
    else
        return MB.algebra_element(res)
    end
end
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

function _multi_basis_type(::Type{BT}, ::Type{T}) where {BT<:SA.SubBasis,T}
    AE = MB.algebra_element_type(
        Vector{T},
        MA.promote_operation(MB.implicit_basis, BT),
    )
    return Vector{MB.SimpleBasis{AE}}
end
function MA.promote_operation(
    ::typeof(SumOfSquares.Certificate.gram_basis),
    C::Type{<:Ideal{SubC}},
) where {SubC}
    return _multi_basis_type(
        MA.promote_operation(SumOfSquares.Certificate.gram_basis, SubC),
        _coeff_type(C),
    )
end
function MA.promote_operation(
    ::typeof(SumOfSquares.Certificate.zero_basis),
    ::Type{<:Ideal{C}},
    ::Type{B},
    ::Type{D},
    ::Type{G},
    ::Type{W},
) where {C,B,D,G,W}
    inner = MA.promote_operation(
        SumOfSquares.Certificate.zero_basis,
        C,
        B,
        D,
        G,
        W,
    )
    return _invariant_basis_type(inner)
end

function _invariant_basis_type(::Type{EB}) where {EB<:SA.ExplicitBasis}
    T = eltype(EB)
    IB = MA.promote_operation(MB.implicit_basis, EB)
    return InvariantBasis{T,Int,IB,EB,SparseArrays.SparseVector{Float64,Int}}
end
_invariant_basis_type(::Type{T}) where {T} = T  # passthrough for non-ExplicitBasis (e.g. QuotientBasis)
SumOfSquares.Certificate.zero_basis(::Ideal) = MB.Monomial
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
function SumOfSquares.Certificate.zero_basis(
    certificate::Ideal,
    basis,
    domain,
    gram_bases,
    weights,
)
    inner = SumOfSquares.Certificate.zero_basis(
        certificate.certificate,
        basis,
        domain,
        gram_bases,
        weights,
    )
    return _invariant_basis(certificate.pattern, inner)
end

function _invariant_basis(pattern::Pattern, monomial_basis::SA.ExplicitBasis)
    tbl = SymbolicWedderburn.Characters.CharacterTable(
        Rational{Int},
        pattern.group,
    )
    invs = SymbolicWedderburn.invariant_vectors(
        tbl,
        pattern.action,
        monomial_basis,
    )
    # Coerce to a uniform `SparseVector{Float64,Int}` so that the bridge's
    # concrete parametric type can be predicted by `promote_operation`.
    invs_uniform = [convert(SparseArrays.SparseVector{Float64,Int}, iv) for iv in invs]
    return InvariantBasis(monomial_basis, invs_uniform)
end
# Fallback for non-SubBasis inner zero_basis (e.g. QuotientBasis): leave it.
_invariant_basis(::Pattern, basis) = basis
function MA.promote_operation(
    ::typeof(SumOfSquares.Certificate.zero_basis),
    ::Type{Ideal{S,C}},
    ::Type{B},
    ::Type{D},
    ::Type{G},
    ::Type{W},
) where {S,C,B,D,G,W}
    return MA.promote_operation(
        SumOfSquares.Certificate.zero_basis,
        C,
        B,
        D,
        G,
        W,
    )
end

function matrix_reps(pattern, R, basis, ::Type{T}, form) where {T}
    monos = MB.keys_as_monomials(basis)
    polys = R * monos
    return map(SymbolicWedderburn.gens(pattern.group)) do g
        S = Matrix{T}(undef, length(polys), length(polys))
        for i in eachindex(polys)
            p = polys[i]
            q = SymbolicWedderburn.action(pattern.action, g, p)
            coefs = MP.coefficients(q, MB.keys_as_monomials(basis))
            S[:, i] = _linsolve(R, coefs, form)
        end
        return S
    end
end

function _coeff_type(C::Type{<:Ideal})
    return SumOfSquares._complex(Float64, SumOfSquares.matrix_cone_type(C))
end

function MB.implicit_basis(certificate::Ideal)
    return MB.implicit_basis(certificate.certificate)
end

function SumOfSquares.Certificate.gram_basis(cert::Ideal, poly)
    basis = SumOfSquares.Certificate.gram_basis(cert.certificate, poly)
    return _gram_basis(cert.pattern, basis, _coeff_type(typeof(cert)))
end

function SumOfSquares.Certificate.gram_weights(
    cert::Ideal,
    gram_basis,
    poly,
    ::Type{T},
) where {T}
    # `gram_basis` is the Pattern result: `Vector{SimpleBasis}` of length |χ|.
    # We multiply each χ-block's contribution by `degree(χ)` so that the
    # invariant-vector-projected SDP constraint matches the polynomial.
    inner_basis = SumOfSquares.Certificate.gram_basis(cert.certificate, poly)
    degrees = gram_character_degrees(cert.pattern, inner_basis, _coeff_type(typeof(cert)))
    # Build per-χ weights via `constant_algebra_element(_, T(d))` so that the
    # SparseCoefficients storage matches what the bridge's `W` parameter
    # expects (Tuple-backed, like the default `MB.constant_algebra_element(_, T)`).
    return [MB.constant_algebra_element(SA.basis(poly), T(d)) for d in degrees]
end

function _fixed_basis(F, basis)
    return MB.SimpleBasis(
        map(eachrow(F)) do row
            ae = MB.implicit(MB.algebra_element(collect(row), basis))
            # Copy coefficients to avoid sharing the basis.keys vector,
            # which would be mutated by `canonical` (called during hash/==)
            return SA.AlgebraElement(copy(SA.coeffs(ae)), Base.parent(ae))
        end,
    )
end

function _gram_basis(pattern::Pattern, basis, ::Type{T}) where {T}
    # `semisimple=false` returns one simple `DirectSummand` per character χ.
    # The numerical block-diagonalization that previously lived here is now
    # inside `SymbolicWedderburn` (numerical_simple.jl) and runs automatically
    # if the symbolic minimal projection cannot reduce to a simple summand.
    summands = SymbolicWedderburn.symmetry_adapted_basis(
        T,
        pattern.group,
        pattern.action,
        basis,
    )
    return map(summands) do summand
        R = convert(Matrix{T}, SymbolicWedderburn.image_basis(summand))
        return _fixed_basis(R, basis)
    end
end

"""
    gram_character_degrees(pattern, basis, ::Type{T})

Return a `Vector{Int}` of degrees of the irreducible characters, parallel to
[`_gram_basis`](@ref). Used to scale each χ-block by `d_χ` so that the
invariant-vector-projected SDP constraint is satisfied (see
`sos_problem.jl` in `SymbolicWedderburn/examples` for the analogous trick).
"""
function gram_character_degrees(pattern::Pattern, basis, ::Type{T}) where {T}
    summands = SymbolicWedderburn.symmetry_adapted_basis(
        T,
        pattern.group,
        pattern.action,
        basis,
    )
    return [SymbolicWedderburn.degree(s) for s in summands]
end
