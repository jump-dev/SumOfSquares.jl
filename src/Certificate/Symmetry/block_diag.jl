import DataStructures

"""
    _reorder!(F::LinearAlgebra.Schur{T}) where {T}

Given a (quasi) upper triangular matrix `S` returns an othogonal
matrix `P` such that `P' * S * P` is still quasi upper triangular
but has its eigenvalues in increasing order.

The matrix `P` is a product of permutation matrices
(to reorder the eigenvalues) and
Givens rotation (to keep it upper triangular after reordering).

By (quasi), we mean that if `S` is a `Matrix{<:Real}`,
then there may be nonzero entries in `S[i+1,i]` representing
complex conjugates.
In that case, the complex conjugate are permuted together.
If `S` is a `Matrix{<:Complex}`, then `S` is triangular.
"""
function _reorder!(F::LinearAlgebra.Schur{T}) where {T}
    n = length(F.values)
    # Bubble sort
    sorted = false
    while !sorted
        prev_i = nothing
        sorted = true
        i = 1
        while i <= n
            S = F.Schur
            if (T <: Real) && i < n && !iszero(S[i+1, i])
                # complex pair
                next_i = i + 2
            else
                next_i = i + 1
            end
            if !isnothing(prev_i) &&
               (real(S[i, i]), imag(S[i, i])) <
               (real(S[prev_i, prev_i]), imag(S[prev_i, prev_i]))
                select = trues(n)
                select[prev_i:(i-1)] .= false
                select[next_i:end] .= false
                LinearAlgebra.ordschur!(F, select)
                sorted = false
            end
            prev_i = i
            i = next_i
        end
    end
end

# We can multiply by `Diagonal(d)` if `d[i] * conj(d[i]) = 1`.
# So in the real case, `d = ±1` but in the complex case, we have more freedom.
function _sign_diag(
    A::AbstractMatrix{T},
    B::AbstractMatrix{T};
    tol = Base.rtoldefault(real(T)),
) where {T}
    n = LinearAlgebra.checksquare(A)
    d = ones(T, n)
    for j in 2:n
        if T <: Real
            minus = zero(real(T))
            not_minus = zero(real(T))
            for i in 1:(j-1)
                a = A[i, j]
                b = B[i, j]
                minus = max(minus, abs(a + b))
                not_minus = max(not_minus, abs(a - b))
            end
            if minus < not_minus
                d[j] = -one(T)
                B[:, j] = -B[:, j]
                B[j, :] = -B[j, :]
            end
        else
            i = argmax(abs.(B[1:(j-1), j]))
            if abs(B[i, j]) <= tol
                continue
            end
            rot = B[i, j] / A[i, j]
            # It should be unitary but there might be small numerical errors
            # so let's normalize
            rot /= abs(rot)
            rot = conj(rot)
            d[j] = rot
            B[:, j] *= rot
            B[j, :] *= rot
        end
    end
    return d
end

# Try to eliminate rounding error if it's easy
# Should be the case for many groups
function _try_integer!(A::Matrix)
    if all(a -> round(a) ≈ a, A)
        return round.(A)
    else
        return A
    end
end

"""
    _rotate_complex(A::AbstractMatrix{T}, B::AbstractMatrix{T}; tol = Base.rtoldefault(real(T))) where {T}

Given (quasi) upper triangular matrix `A` and `B` that have the eigenvalues in
the same order except the complex pairs which may need to be (signed) permuted,
returns an othogonal matrix `P` such that `P' * A * P` and `B` have matching
low triangular part.
The upper triangular part will be dealt with by `_sign_diag`.

By (quasi), we mean that if `S` is a `Matrix{<:Real}`,
then there may be nonzero entries in `S[i+1,i]` representing
complex conjugates.
If `S` is a `Matrix{<:Complex}`, then `S` is upper triangular so there is
nothing to do.
"""
function _rotate_complex(
    A::AbstractMatrix{T},
    B::AbstractMatrix{T};
    tol = Base.rtoldefault(real(T)),
) where {T}
    n = LinearAlgebra.checksquare(A)
    I = collect(1:n)
    J = copy(I)
    V = ones(T, n)
    pair = false
    for i in 1:n
        if pair || i == n
            continue
        end
        pair = abs(A[i+1, i]) > tol
        if pair
            a = (A[i+1, i], A[i, i+1])
            b = (B[i+1, i], B[i, i+1])
            c = a[2:-1:1]
            if LinearAlgebra.norm(abs.(a) .- abs.(b)) >
               LinearAlgebra.norm(abs.(c) .- abs.(b))
                a = c
                J[i] = i + 1
                J[i+1] = i
            end
            c = (-).(a)
            if LinearAlgebra.norm(a .- b) > LinearAlgebra.norm(c .- b)
                V[i+1] = -V[i]
            end
        end
    end
    return SparseArrays.sparse(I, J, V, n, n)
end

"""
    orthogonal_transformation_to(A, B)

Return an orthogonal transformation `U` such that
`A = U' * B * U`

Given Schur decompositions
`A = Z_A * S_A * Z_A'`
`B = Z_B * S_B * Z_B'`
Since `P' * S_A * P = S_B`, we have
`A = Z_A * P * Z_B' * B * Z_B * P' * Z_A'`
"""
function orthogonal_transformation_to(A, B)
    n = LinearAlgebra.checksquare(A)
    As = LinearAlgebra.schur(A)
    _reorder!(As)
    T_A = As.Schur
    Z_A = As.vectors
    Bs = LinearAlgebra.schur(B)
    _reorder!(Bs)
    T_B = Bs.Schur
    Z_B = Bs.vectors
    P = _rotate_complex(T_A, T_B)
    T_A = P' * T_A * P
    d = _sign_diag(T_A, T_B)
    return _try_integer!(Z_B * LinearAlgebra.Diagonal(d) * P' * Z_A')
end

function ordered_block_diag(As, d)
    U = block_diag(As, d)
    U === nothing && return nothing
    iU = U'
    @assert iU ≈ inv(U)
    Bs = [iU * A * U for A in As]
    @assert all(Bs) do B
        return isblockdim(B, d)
    end
    refs = [B[1:d, 1:d] for B in Bs]
    offset = d
    for offset in d:d:(size(U, 1)-d)
        I = offset .+ (1:d)
        Cs = [B[I, I] for B in Bs]
        λ = rand(length(Bs))
        # We want to find a transformation such that
        # the blocks `Cs` are equal to the blocks `refs`
        # With probability one, making a random combination match
        # should work, this trick is similar to [CGT97].
        #
        # [CGT97] Corless, R. M.; Gianni, P. M. & Trager, B. M.
        # A reordered Schur factorization method for zero-dimensional polynomial systems with multiple roots Proceedings of the 1997 international symposium on Symbolic and algebraic computation,
        # 1997, 133-140
        R = sum(λ .* refs)
        C = sum(λ .* Cs)
        V = orthogonal_transformation_to(R, C)
        @assert R ≈ V' * C * V
        for i in eachindex(refs)
            @assert refs[i] ≈ V' * Cs[i] * V
        end
        U[:, I] = U[:, I] * V
        offset += d
    end
    @assert ordered_block_check(U, As, d)
    return U
end

function ordered_block_check(U, As, d)
    iU = U'
    if !(iU ≈ inv(U))
        return false
    end
    return all(As) do A
        return is_ordered_blockdim(iU * A * U, d)
    end
end

function block_diag(As, d)
    for A in As
        #T = LinearAlgebra.eigen(A).vectors
        Z = LinearAlgebra.schur(A).vectors
        iZ = Z'
        @assert iZ ≈ inv(Z)
        n = LinearAlgebra.checksquare(A)
        union_find = DataStructures.IntDisjointSets(n)
        Bs = [iZ * A * Z for A in As]
        for B in Bs
            merge_sparsity!(union_find, B)
        end
        blocks = Dict{Int,Vector{Int}}()
        for i in 1:n
            r = DataStructures.find_root!(union_find, i)
            blocks[r] = push!(get(blocks, r, Int[]), i)
        end
        if length(blocks) > 1
            U = similar(Z)
            offset = 0
            for v in values(blocks)
                @assert iszero(length(v) % d)
                if length(v) == d
                    V = Z[:, v]
                else
                    Cs = [B[v, v] for B in Bs]
                    V = block_diag(Cs, d)
                    V === nothing && break
                    V *= transpose(Z[:, v])
                end
                U[:, offset.+eachindex(v)] = V
                offset += length(v)
            end
            if offset == n
                return U
            end
        end
    end
end

function merge_sparsity!(
    union_find::DataStructures.IntDisjointSets,
    A,
    tol = 1e-8,
)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(A[I]) > tol
            union!(union_find, i, j)
        end
    end
end

function is_ordered_blockdim(A, d, tol = 1e-8)
    B = A[1:d, 1:d]
    return isblockdim(A, d, tol) && all(d:d:(size(A, 1)-d)) do offset
        I = offset .+ (1:d)
        return isapprox(B, A[I, I], rtol = tol)
    end
end
function isblockdim(A, d, tol = 1e-8)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(i - j) >= d && abs(A[I]) > tol
            return false
        end
    end
    return true
end
