"""
    struct InvariantBasis{T,I,IB,EB,IV} <: SA.ExplicitBasis{T,I}
        implicit_basis::IB     # parent implicit basis (e.g. MB.FullBasis{Monomial})
        monomial_basis::EB     # SubBasis over which invariant_vectors index
        invariant_vectors::Vector{IV}
    end

An explicit basis where each "element" represents a G-invariant orbit of
monomials, encoded as a sparse vector over `monomial_basis`. Coefficients of a
polynomial in this basis are the projections of the monomial coefficients onto
each invariant vector.

Used by `SumOfSquares.Certificate.Symmetry` so that the SDP constraint side
emits one scalar equality per invariant orbit rather than one per monomial,
which is what allows the gram side to use a single simple basis per irreducible
character without breaking polynomial-identity equality.
"""
struct InvariantBasis{T,I,IB<:SA.ImplicitBasis,EB<:SA.ExplicitBasis,IV<:SparseArrays.AbstractSparseVector} <:
       SA.ExplicitBasis{T,I}
    implicit_basis::IB
    monomial_basis::EB
    invariant_vectors::Vector{IV}
end

function InvariantBasis(
    monomial_basis::SA.ExplicitBasis,
    invariant_vectors::Vector{<:SparseArrays.AbstractSparseVector},
)
    IB = parent(monomial_basis)
    T = eltype(IB)
    return InvariantBasis{T,Int,typeof(IB),typeof(monomial_basis),eltype(invariant_vectors)}(
        IB,
        monomial_basis,
        invariant_vectors,
    )
end

Base.length(b::InvariantBasis) = length(b.invariant_vectors)
Base.parent(b::InvariantBasis) = b.implicit_basis
MB.implicit_basis(b::InvariantBasis) = b.implicit_basis

function Base.:(==)(a::InvariantBasis, b::InvariantBasis)
    return a.implicit_basis == b.implicit_basis &&
           a.monomial_basis == b.monomial_basis &&
           a.invariant_vectors == b.invariant_vectors
end

# Override `coeffs(cfs, source, target::InvariantBasis)` so that translating
# coefficients into an invariant basis performs the invariant-vector projection.
function SA.coeffs(cfs, source::SA.AbstractBasis, target::InvariantBasis)
    res = SA.zero_coeffs(SA.value_type(cfs), target)
    return SA.coeffs!(res, cfs, source, target)
end

function SA.coeffs!(res, cfs, source::SA.AbstractBasis, target::InvariantBasis)
    MA.operate!(zero, res)
    mb = target.monomial_basis
    for (k, v) in SA.nonzero_pairs(cfs)
        mono = source[k]
        # `get(::SubBasis, key, ::Nothing)` returns the integer position or `nothing`.
        m = get(mb, mono, nothing)
        isnothing(m) && continue
        for ki in eachindex(target.invariant_vectors)
            iv_at_m = target.invariant_vectors[ki][m]
            iszero(iv_at_m) && continue
            res[ki] = MA.operate!!(MA.add_mul, res[ki], iv_at_m, v)
        end
    end
    return res
end
