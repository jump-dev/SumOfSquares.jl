mutable struct XORSpace{T<:Integer}
    dimension::Int
    basis::Vector{T}
    function XORSpace{T}(n) where {T}
        return new(0, zeros(T, n))
    end
end
function appropriate_type(n)
    if n < 8sizeof(Int)
        return Int
    elseif n < 64
        return Int64
    elseif n < 128
        return Int128
    else
        return BigInt
    end
end
function XORSpace(n)
    return XORSpace{appropriate_type(n)}(n)
end
function Base.push!(xs::XORSpace{T}, x::T) where {T}
    for i in eachindex(xs.basis)
        if !iszero(x & (one(T) << (i - 1)))
            if iszero(xs.basis[i])
                xs.dimension += 1
                xs.basis[i] = x
                break
            else
                x = xor(x, xs.basis[i])
            end
        end
    end
    return xs
end
function orthogonal_complement(xs::XORSpace{T}) where {T}
    r = Vector{T}(undef, length(xs.basis) - xs.dimension)
    k = 0
    for i in eachindex(xs.basis)
        if iszero(xs.basis[i])
            e = one(T) << (i - 1)
            x = e
            for j in 1:(i-1)
                if !iszero(xs.basis[j] & e)
                    x |= (one(T) << (j - 1))
                end
            end
            k += 1
            r[k] = x
        end
    end
    @assert k == length(r)
    return r
end
function xor_complement(x, n, ::Type{T} = appropriate_type(n)) where {T}
    xs = XORSpace{T}(n)
    for xi in x
        push!(xs, xi)
    end
    return orthogonal_complement(xs)
end
