# Inspired from SOSTools
export randpsd, randsos

function randpsd(n; r=n, eps=0.1)
    Q = randn(n,n)
    d = zeros(Float64, n)
    d[1:r] = eps .+ abs.(randn(r))
    return Q' * Diagonal(d) * Q
end

function _randsos(X::AbstractVector{<:MP.AbstractMonomial}; r=-1, monotype=:Classic, eps=0.1)
    if monotype == :Classic
        x = Certificate.monomials_half_newton_polytope(X, tuple())
    elseif monotype == :Gram
        x = X
    else
        throw(ArgumentError("Monotype $monotype not known"))
    end
    n = length(x)
    if r < 0
        r = n
    end
    return GramMatrix(randpsd(n, r=r, eps=eps), x)
end

randsos(X::AbstractVector; kws...) = _randsos(MP.monovec(X); kws...)
