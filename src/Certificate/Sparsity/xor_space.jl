"""
    mutable struct XORSpace{T<:Integer}
        dimension::Int
        basis::Vector{T}
    end

Basis for a linear subspace of the Hamming space (i.e. set of binary string
`{0, 1}^n` of length `n`) represented in the bits of an integer of type `T`.
This is used to implement `Certificate.Sparsity.SignSymmetry`.

Consider the scalar product `⟨a, b⟩` which returns the
the xor of the bits of `a & b`.
It is a scalar product since `⟨a, b⟩ = ⟨b, a⟩` and
`⟨a, xor(b, c)⟩ = xor(⟨a, b⟩, ⟨a, c⟩)`.

We have two options here to compute the orthogonal space.

The first one is to build an orthogonal basis
with some kind of Gram-Schmidt process and then to obtain the orthogonal
space by removing the projection from each vector of the basis.

The second option, which is the one we use here is to compute the row echelon
form and then read out the orthogonal subspace directly from it.
For instance, if the row echelon form is
```
1 a 0 c e
  b 1 d f
```
then the orthogonal basis is
```
a 1 b 0 0
c 0 d 1 0
e 0 f 0 1
```

The `basis` vector has `dimension` nonzero elements. Any element added with
`push!` can be obtained as the `xor` of some of the elements of `basis`.
Moreover, the `basis` is kept in row echelon form in the sense that the first,
second, ..., `i - 1`th bits of `basis[i]` are zero and `basis[i]` is zero if its
`i`th bit is not one.
"""
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
b(x) = bitstring(x)[end-5:end]
e_i(T, i) = (one(T) << (i - 1))
has_bit(x, i) = !iszero(x & e_i(typeof(x), i))
function Base.push!(xs::XORSpace{T}, x::T) where {T}
    I = eachindex(xs.basis)
    for i in I
        if !iszero(xs.basis) && has_bit(x, i)
            x = xor(x, xs.basis[i])
        end
    end
    i = findfirst(Base.Fix1(has_bit, x), I)
    if isnothing(i)
        return
    end
    xs.dimension += 1
    # Clear `i`th column elements of the basis already added
    for j in I
        if has_bit(xs.basis[j], i)
            xs.basis[j] = xor(xs.basis[j], x)
        end
    end
    # `x` will be the row used as pivot for the `i`th column
    xs.basis[i] = x
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
