function is_orthogonal(S, tol=1e-6)
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
    return LinearAlgebra.dot(v, w) / LinearAlgebra.dot(w, w)
end
function _proj(A, i, v)
    w = A[i, :]
    return w * _projcoef(A, i, v, w)
end

function orthogonalize(A)
    display(A)
    for i in 2:size(A, 1)
        row = A[i, :]
        for j in 1:(i - 1)
            row = row - _proj(A, j, row)
        end
        A[i, :] = row
    end
    return A
end

struct _OrthogonalMatrix end
struct _RowEchelonMatrix end

function __linsolve(A::Matrix{T}, b::Vector{T}, ::_OrthogonalMatrix) where {T}
    return map(1:size(A, 1)) do i
        return _projcoef(A, i, b)
    end
end

function __linsolve(A::Matrix{T}, b::Vector{T}, ::_RowEchelonMatrix) where {T}
    j = 0
    return map(1:size(A, 1)) do i
        while j < size(A, 2)
            j += 1
            if isone(A[i, j]) && all(k -> k == i || iszero(A[k, j]), 1:size(A, 1))
                return b[j]
            end
        end
        error("Not in row_echelon_form, cannot find for `$i`th entry.")
    end
end

function _linsolve(A::Matrix{T}, b::Vector{T}, form) where {T}
    x = __linsolve(A, b, form)
    #@assert transpose(A) * x == b
    return x
end
