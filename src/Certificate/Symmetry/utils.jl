function is_orthogonal(S, tol = 1e-6)
    n = LinearAlgebra.checksquare(S)
    for i in 1:n
        for j in 1:n
            if i != j
                s = LinearAlgebra.dot(S[:, i], S[:, j])
                if abs(s) > tol
                    return false
                end
            end
        end
    end
    return true
end

function _projcoef(A, i, v, w = A[i, :])
    # `LinearAlgebra.dot(v, w)` would not work as it is computing the adjoint
    # of `v`.
    # If `v = α * w + ...` then
    # `dot(w, v) = conj(w) (α * w + ...) = α`
    # However,
    # `dot(v, w) = conj(α * w + ...) w = conj(α)`
    return LinearAlgebra.dot(w, v) / LinearAlgebra.dot(w, w)
end

struct _OrthogonalMatrix end
struct _RowEchelonMatrix end

function __linsolve(
    A::AbstractMatrix{T},
    b::Vector{T},
    ::_OrthogonalMatrix,
) where {T}
    return map(1:size(A, 1)) do i
        return _projcoef(A, i, b)
    end
end

function __linsolve(A::Matrix{T}, b::Vector{T}, ::_RowEchelonMatrix) where {T}
    j = 0
    return map(1:size(A, 1)) do i
        while j < size(A, 2)
            j += 1
            if isone(A[i, j]) &&
               all(k -> k == i || iszero(A[k, j]), 1:size(A, 1))
                return b[j]
            end
        end
        return error("Not in row_echelon_form, cannot find for `$i`th entry.")
    end
end

using SparseArrays
function _linsolve(A::AbstractMatrix{T}, b::Vector{T}, form) where {T}
    x = __linsolve(A, b, form)
    @assert transpose(A) * x ≈ b
    return x
end
