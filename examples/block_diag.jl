using LinearAlgebra
using Combinatorics, DataStructures

function row_echelon_linsolve(A::Matrix{T}, b::Vector{T}) where {T}
    j = 0
    x = map(1:size(A, 1)) do i
        while j < size(A, 2)
            j += 1
            if isone(A[i, j]) && all(k -> k == i || iszero(A[k, j]), 1:size(A, 1))
                return b[j]
            end
        end
        error("Not in row_echelon_form, cannot find for `$i`th entry.")
    end
    @assert x'A == b'
    return x
end

function ordered_block_diag(As, d)
    U = block_diag(As, d)
    U === nothing && return nothing
    iU = inv(U)
    Bs = [U * A * inv(U) for A in As]
    @assert all(Bs) do B
        isblockdim(B, d)
    end
    refs = [B[1:d, 1:d] for B in Bs]
    offset = d
    for offset in d:d:(size(U, 1) - d)
        I = offset .+ (1:d)
        Cs = [B[I, I] for B in Bs]
        σ_ok = nothing
        for σ in permutations(1:d)
            if all(zip(refs, Cs)) do refC
                ref, C = refC
                isapprox(ref[σ, σ], C, rtol=1e-8)
            end
                σ_ok = σ
                break
            end
        end
        if σ_ok === nothing
            error("No permutation can make $refs and $Cs match")
        end
        U[:, I] = U[:, I[σ_ok]]
        offset += d
    end
    iU = inv(U)
    @assert all(As) do A
        is_ordered_blockdim(U * A * iU, d)
    end
    return U
end

function block_diag(As, d)
    for A in As
        #T = eigen(A).vectors
        T = schur(A).vectors
        iT = inv(T)
        n = LinearAlgebra.checksquare(A)
        union_find = IntDisjointSets(n)
        Bs = [iT * A * T for A in As]
        for B in Bs
            merge_sparsity!(union_find, B)
        end
        blocks = Dict{Int,Vector{Int}}()
        for i in 1:n
            r = find_root!(union_find, i)
            blocks[r] = push!(get(blocks, r, Int[]), i)
        end
        if length(blocks) > 1
            U = similar(T)
            offset = 0
            for v in values(blocks)
                @assert iszero(length(v) % d)
                if length(v) == d
                    V = T[:, v]
                else
                    Cs = [B[v, v] for B in Bs]
                    V = block_diag(Cs, d)
                    V === nothing && break
                    V *= T[:, v]'
                end
                U[:, offset .+ eachindex(v)] = V
                offset += length(v)
            end
            if offset == n
                return U
            end
        end
    end
end

function merge_sparsity!(union_find::IntDisjointSets, A, tol=1e-8)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(A[I]) > tol
            union!(union_find, i, j)
        end
    end
end

function is_ordered_blockdim(A, d, tol=1e-8)
    B = A[1:d, 1:d]
    return isblockdim(A, d, tol) &&
        all(d:d:(size(A, 1) - d)) do offset
            I = offset .+ (1:d)
            isapprox(B, A[I, I], rtol=tol)
        end
end
function isblockdim(A, d, tol=1e-8)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(i - j) >= d && abs(A[I]) > tol
            return false
        end
    end
    return true
end
