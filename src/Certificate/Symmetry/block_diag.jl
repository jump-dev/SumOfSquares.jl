import Combinatorics, DataStructures

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
        σ_ok = nothing
        for σ in Combinatorics.permutations(1:d)
            if all(zip(refs, Cs)) do refC
                ref, C = refC
                return isapprox(ref[σ, σ], C, rtol = 1e-8)
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
    ordered_block_check(U, As, d)
    return U
end

function ordered_block_check(U, As, d)
    iU = U'
    @assert iU ≈ inv(U)
    @assert all(As) do A
        return is_ordered_blockdim(iU * A * U, d)
    end
end

function block_diag(As, d)
    for A in As
        #T = LinearAlgebra.eigen(A).vectors
        T = LinearAlgebra.schur(A).vectors
        iT = T'
        @assert iT ≈ inv(T)
        n = LinearAlgebra.checksquare(A)
        union_find = DataStructures.IntDisjointSets(n)
        Bs = [iT * A * T for A in As]
        for B in Bs
            merge_sparsity!(union_find, B)
        end
        blocks = Dict{Int,Vector{Int}}()
        for i in 1:n
            r = DataStructures.find_root!(union_find, i)
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
                    V *= transpose(T[:, v])
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
