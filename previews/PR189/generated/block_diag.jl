using LinearAlgebra
using DataStructures

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

function block_diag(As, m)
    for A in As
        #T = eigen(A).vectors
        T = schur(A).vectors
        #@show D
        iT = inv(T)
        n = LinearAlgebra.checksquare(A)
        union_find = IntDisjointSets(n)
        Bs = [iT * A * T for A in As]
        for B in Bs
            merge_sparsity!(union_find, B)
            #@show union_find
        end
        d = Dict{Int,Vector{Int}}()
        for i in 1:n
            r = find_root!(union_find, i)
            d[r] = push!(get(d, r, Int[]), i)
        end
        #@show d
        if length(d) > 1
            U = similar(T)
            offset = 0
            for v in values(d)
                @assert iszero(length(v) % m)
                if length(v) == m
                    V = T[:, v]
                else
                    Cs = [B[v, v] for B in Bs]
                    V = block_diag(Cs, m)
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

function isblockdim(A, m, tol=1e-8)
    for I in CartesianIndices(A)
        i, j = I.I
        if abs(i - j) >= m && abs(A[I]) > tol
            return false
        end
    end
    return true
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

