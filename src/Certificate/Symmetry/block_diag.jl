import DataStructures

"""
    _permutation_quasi_upper_triangular(S)

Given a (quasi) upper triangular matrix `S`
returns the permutation `P` so that
`P' * S * P` has its eigenvalues in increasing order.
"""
function _permutation_quasi_upper_triangular(S::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(S)
    # Bubble sort
    sorted = false
    P = SparseArrays.sparse(one(T) * LinearAlgebra.I, n, n)
    function permute!(i, j)
        swap = sparse([i, j], [j, i], ones(T, 2), n, n)
        S = swap' * S * swap
        return P *= swap
    end
    while !sorted
        prev_i = nothing
        sorted = true
        i = 1
        while i <= n
            if (i < n && !iszero(S[i+1, i]))
                #if S[i+1, i] < S[i, i+1]
                #    permute!(i, i + 1)
                #end
                if !isnothing(prev_i) && S[i, i] < S[prev_i, prev_i]
                    if i - prev_i == 2
                        permute!(prev_i, i)
                        permute!(prev_i + 1, i + 1)
                    else
                        permute!(prev_i, i)
                        permute!(i, i + 1)
                    end
                    sorted = false
                end
                # complex
                prev_i = i
                i += 2
            else
                if !isnothing(prev_i) && S[i, i] < S[prev_i, prev_i]
                    permute!(prev_i, i)
                    if i - prev_i == 2
                        permute!(i - 1, i)
                    end
                    sorted = false
                end
                prev_i = i
                i += 1
            end
        end
    end
    return P
end

# `A` and `B` may be not upper triangular because of the permutations
function _sign_diag(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)
    d = ones(T, n)
    for i in 1:n
        minus = zero(T)
        not_minus = zero(T)
        for j in 1:(i-1)
            for (I, J) in [(i, j), (j, i)]
                a = A[I, J]
                b = B[I, J]
                minus = max(minus, abs(a + b))
                not_minus = max(not_minus, abs(a - b))
            end
        end
        if minus < not_minus
            d[i] = -one(T)
            B[:, i] = -B[:, i]
            B[i, :] = -B[i, :]
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
    orthogonal_transformation_to(A, B)

Return an orthogonal transformation `U` such that
`A = U' * B * U`

Given Schur decompositions
`A = Z_A * S_A * Z_A'`
`B = Z_B * S_B * Z_B'`
We further decompose the triangular matrices `S_A`, `S_B`
to order the eigenvalues:
`S_A = P_A * T_A * P_A'`
`S_B = P_B * T_B * P_B'`
"""
function orthogonal_transformation_to(A, B)
    As = LinearAlgebra.schur(A)
    T_A = As.Schur
    Z_A = As.vectors
    P_A = _permutation_quasi_upper_triangular(T_A)
    Bs = LinearAlgebra.schur(B)
    T_B = Bs.Schur
    Z_B = Bs.vectors
    P_B = _permutation_quasi_upper_triangular(T_B)
    d = _sign_diag(P_A' * T_A * P_A, P_B' * T_B * P_B)
    return _try_integer!(Z_B * P_B * LinearAlgebra.Diagonal(d) * P_A' * Z_A')
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
