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
    buckets = Dict{U, Vector{eltype(monos)}}()
    for mono in monos
        exp = binary_exponent(MP.exponents(mono), T)
        i = 1 + binary_exponent([bin_dot(ri, exp) for ri in r], U)
        if !haskey(buckets, i)
            buckets[i] = eltype(monos)[]
        end
        push!(buckets[i], mono)
    end
    return collect(values(buckets))
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
