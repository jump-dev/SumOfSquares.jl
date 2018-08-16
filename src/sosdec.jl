export SOSDecomposition

struct SOSDecomposition{T, PT <: AbstractPolynomialLike{T}} <: AbstractPolynomialLike{T} # If promote_op((x, y) -> x * y + x * y, T, T) != T then it might not be true
    ps::Vector{PT}
    function SOSDecomposition{T, PT}(ps::Vector{PT}) where {T, PT}
        new(ps)
    end
end
SOSDecomposition(ps::Vector{PT}) where {T, PT <: AbstractPolynomialLike{T}} = SOSDecomposition{T, PT}(ps)
MP.polynomialtype(::Type{SOSDecomposition{T, PT}}) where {T, PT} = polynomialtype(PT)
#function SOSDecomposition(ps::Vector)
#    T = reduce(promote_type, Int, map(eltype, ps))
#    SOSDecomposition{T}(ps)
#end

function MatPolynomial(p::SOSDecomposition{T}) where {T}
    X = mergemonovec(map(monomials, p))
    m = length(p)
    n = length(X)
    Q = zeros(T, m, n)
    for i in 1:m
        j = 1
        for t in terms(p[i])
            while X[j] != monomial(t)
                j += 1
            end
            Q[i, j] = coefficient(t)
            j += 1
        end
    end
    MatPolynomial(Q' * Q, X)
end

function SOSDecomposition(p::MatPolynomial)
    n = length(p.x)
    # TODO LDL^T factorization for SDP is missing in Julia
    # it would be nice to have though
    A = getmat(p)
    @static if VERSION >= v"0.7-"
        Q = cholesky(Matrix(A)).U
    else
        Q = chol(A)
    end
    m = size(Q, 1)
    ps = [polynomial(Q[i,:], p.x) for i in 1:m]
    SOSDecomposition(ps)
end
# Without LDL^T, we need to do float(T)
SOSDecomposition(p::MatPolynomial{C, T}) where {C, T} = SOSDecomposition{C, float(T)}(p)

Base.length(p::SOSDecomposition) = length(p.ps)
Base.isempty(p::SOSDecomposition) = isempty(p.ps)
Base.start(p::SOSDecomposition) = start(p.ps)
Base.done(p::SOSDecomposition, state) = done(p.ps, state)
Base.next(p::SOSDecomposition, state) = next(p.ps, state)
Base.getindex(p::SOSDecomposition, i::Int) = p.ps[i]

(p::MatPolynomial)(s::MP.AbstractSubstitution...) = polynomial(p)(s...)

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
