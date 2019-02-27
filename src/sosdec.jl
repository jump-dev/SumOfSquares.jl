export SOSDecomposition

struct SOSDecomposition{T, PT <: MP.APL{T}} <: MP.APL{T} # If promote_op((x, y) -> x * y + x * y, T, T) != T then it might not be true
    ps::Vector{PT}
    function SOSDecomposition{T, PT}(ps::Vector{PT}) where {T, PT}
        new(ps)
    end
end
SOSDecomposition(ps::Vector{PT}) where {T, PT <: MP.APL{T}} = SOSDecomposition{T, PT}(ps)
MP.polynomialtype(::Type{SOSDecomposition{T, PT}}) where {T, PT} = MP.polynomialtype(PT)
#function SOSDecomposition(ps::Vector)
#    T = reduce(promote_type, Int, map(eltype, ps))
#    SOSDecomposition{T}(ps)
#end

function GramMatrix(p::SOSDecomposition{T}) where {T}
    X = MP.mergemonovec(map(MP.monomials, p))
    m = length(p)
    n = length(X)
    Q = zeros(T, m, n)
    for i in 1:m
        j = 1
        for t in MP.terms(p[i])
            while X[j] != MP.monomial(t)
                j += 1
            end
            Q[i, j] = MP.coefficient(t)
            j += 1
        end
    end
    GramMatrix(Q' * Q, X)
end

function SOSDecomposition(p::GramMatrix)
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    A = getmat(p)
    Q = cholesky(Matrix(A)).U
    m = size(Q, 1)
    ps = [MP.polynomial(Q[i,:], p.x) for i in 1:m]
    SOSDecomposition(ps)
end
# Without LDL^T, we need to do float(T)
SOSDecomposition(p::GramMatrix{C, T}) where {C, T} = SOSDecomposition{C, float(T)}(p)

Base.length(p::SOSDecomposition) = length(p.ps)
Base.isempty(p::SOSDecomposition) = isempty(p.ps)
Base.iterate(p::SOSDecomposition, args...) = Base.iterate(p.ps, args...)
Base.getindex(p::SOSDecomposition, i::Int) = p.ps[i]

(p::GramMatrix)(s::MP.AbstractSubstitution...) = MP.polynomial(p)(s...)

function Base.show(io::IO, p::SOSDecomposition)
    for (i, q) in enumerate(p)
        print(io, "(")
        print(io, q)
        print(io, ")^2")
        if i != length(p)
            print(io, " + ")
        end
    end
end

function Base.isapprox(p::SOSDecomposition, q::SOSDecomposition; kwargs...)
    m = length(p.ps)
    if length(q.ps) != m
        false
    else
        MultivariateMoments.permcomp((i, j) -> isapprox(p.ps[i], q.ps[j]; kwargs...), m)
    end
end
