"""
    struct SignSymmetry <: Sparsity.Pattern end

Sign symmetry as developed in [Section III.C, L09].
Let `n` be the number of variables.
A *sign-symmetry* is a binary vector `r` of `{0, 1}^n`
such that `dot(r, e)` is even for all exponent `e`.

Let `o(e)` be the binary vector of `{0, 1}^n` for which the `i`th bit is `i`
iff the `i`th entry of `e` is odd.
Let `O` be the set of `o(e)` for exponents of `e`.
The sign symmetry of `r` is then equivalent to its orthogonality
with all elements of `O`.

Since we are only interested in the orthogonal subspace, say `R`, of `O`,
even if `O` is not a linear subspace (i.e., it is not invariant under `xor`),
we compute its span.
We start by creating row echelon form of the span of `O` using
[`XORSpace`](@ref).
We then compute a basis for `R`.
The set `R` of sign symmetries will be the span of this basis.

Let `a`, `b` be exponents of monomials of the gram basis. For any `r in R`,
```
⟨r, o(a + b)⟩ = ⟨r, xor(o(a), o(b)⟩ = xor(⟨r, o(a)⟩, ⟨r, o(b)⟩)
```
For `o(a, b)` to be sign symmetric, this scalar product should be zero for all
sign symmetry `r`.
This is equivalent to saying that `⟨r, o(a)⟩` and `⟨r, o(b)⟩` are equal
for all `r in R`.
In other words, the projection of `o(a)` and `o(b)` to `R` have the same
coordinates in the basis.
If we order the monomials by grouping them by equal coordinates of projection,
we see that the product that are sign symmetric form a block diagonal
structure.
This means that we can group the monomials by these coordinates.

[L09] Löfberg, Johan.
*Pre-and post-processing sum-of-squares programs in practice*.
IEEE transactions on automatic control 54, no. 5 (2009): 1007-1011.
"""
struct SignSymmetry <: Pattern end

function binary_exponent(exponents, ::Type{T}) where {T}
    cur = zero(T)
    for exponent in exponents
        cur <<= 1
        if isodd(exponent)
            cur |= one(T)
        end
    end
    return cur
end
bin_dot(x, y) = isodd(count_ones(x & y))
function buckets_sign_symmetry(monos, r, ::Type{T}, ::Type{U}) where {T,U}
    # We store vector of indices instead of vector of monos in the bucket
    # as `monos` is already sorted and the buckets content will be sorted
    # so it will remain sorted and `DynamicPolynomials` can keep the
    # `MonomialVector` struct in `monos[I]`. With MultivariateBases, the
    # sorted state will be ensured by the basis struct so it will also
    # serve other MultivariatePolynomials implementations
    buckets = Dict{U,Vector{Int}}()
    for (idx, mono) in enumerate(monos)
        exp = binary_exponent(MP.exponents(mono), T)
        i = 1 + binary_exponent([bin_dot(ri, exp) for ri in r], U)
        if !haskey(buckets, i)
            buckets[i] = Int[]
        end
        push!(buckets[i], idx)
    end
    return [monos[I] for I in values(buckets)]
end
function sign_symmetry(
    monos::AbstractVector{<:MP.AbstractMonomial},
    n,
    ::Type{T},
    gram_monos,
) where {T}
    r = xor_complement(
        (binary_exponent(MP.exponents(mono), T) for mono in monos),
        n,
        T,
    )
    return buckets_sign_symmetry(gram_monos, r, T, appropriate_type(length(r)))
end
function sparsity(
    monos::AbstractVector{<:MP.AbstractMonomial},
    ::SignSymmetry,
    gram_monos::AbstractVector{<:MP.AbstractMonomial},
)
    n = MP.nvariables(monos)
    return sign_symmetry(monos, n, appropriate_type(n), gram_monos)
end
