# Inspired from SOSTools
export randpsd, randsos

function randpsd(n; r = n, eps = 0.1)
    Q = randn(n, n)
    d = zeros(Float64, n)
    d[1:r] = eps .+ abs.(randn(r))
    return Q' * Diagonal(d) * Q
end

function _monomials_half_newton_polytope(_monos, filter)
    monos = MP.monomial_vector(_monos)
    basis = MB.FullBasis{MB.Monomial,eltype(monos)}()
    return SumOfSquares.Certificate._half_newton_polytope(
        MB.algebra_element(
            SA.SparseCoefficients(monos, ones(length(monos))),
            basis,
        ),
        MP.variables(monos),
        filter,
    ).monomials
end

function _randsos(
    X::AbstractVector{<:MP.AbstractMonomial};
    r = -1,
    monotype = :Classic,
    eps = 0.1,
)
    if monotype == :Classic
        x = _monomials_half_newton_polytope(
            X,
            Certificate.NewtonDegreeBounds(tuple()),
        )
    elseif monotype == :Gram
        x = X
    else
        throw(ArgumentError("Monotype $monotype not known"))
    end
    n = length(x)
    if r < 0
        r = n
    end
    return GramMatrix(randpsd(n, r = r, eps = eps), x)
end

randsos(X::AbstractVector; kws...) = _randsos(MP.monomial_vector(X); kws...)
