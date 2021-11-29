"""
    struct SignSymmetry <: Sparsity end

Sign symmetry as developed in [Section III.C, L09].

[L09] Lofberg, Johan.
*Pre-and post-processing sum-of-squares programs in practice*.
IEEE transactions on automatic control 54, no. 5 (2009): 1007-1011.
"""
struct SignSymmetry <: Sparsity end

function binary_exponent(exponents, ::Type{T}) where T
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
function buckets_sign_symmetry(monos, r, ::Type{T}, ::Type{U}) where {T, U}
    # We store vector of indices instead of vector of monos in the bucket
    # as `monos` is already sorted and the buckets content will be sorted
    # so it will remain sorted and `DynamicPolynomials` can keep the
    # `MonomialVector` struct in `monos[I]`.
    buckets = Dict{U, Vector{Int}}()
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
function sign_symmetry(monos::AbstractVector{<:MP.AbstractMonomial}, n, ::Type{T}, gram_monos) where T
    r = xor_complement((binary_exponent(MP.exponents(mono), T) for mono in monos), n, T)
    return buckets_sign_symmetry(gram_monos,
                                 r, T, appropriate_type(length(r)))
end
function sparsity(monos::AbstractVector{<:MP.AbstractMonomial}, sp::SignSymmetry,
                  gram_monos::AbstractVector = SumOfSquares.Certificate.monomials_half_newton_polytope(monos, tuple()))
    n = MP.nvariables(monos)
    return sign_symmetry(monos, n, appropriate_type(n), gram_monos)
end
