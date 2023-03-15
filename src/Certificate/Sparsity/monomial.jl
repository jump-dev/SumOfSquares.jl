const CEG = ChordalExtensionGraph

"""
    struct Sparsity.Monomial{C<:CEG.AbstractCompletion} <: Sparsity.Pattern
        completion::C
        k::Int
        use_all_monomials::Bool
    end

Monomial or term sparsity as developed in [WML20a, WML20b].
The `completion` field should be `ClusterCompletion()` [default] for the block-closure or cluster completion [WML20a],
and `ChordalCompletion()` for chordal completion [WML20b].
The integer `k` [default=0] corresponds to `Σ_k` defined in [(3.2), WML20a]
and `k = 0` corresponds to `Σ_*` defined in [(3.3), WML20a].
If `use_all_monomials` is `false` then some monomials of the basis
might be dropped from the basis if not needed.

[WML20a] Wang, Jie, Victor Magron, and Jean-Bernard Lasserre.
*TSSOS: A Moment-SOS hierarchy that exploits term sparsity*.
arXiv preprint arXiv:1912.08899 (2020).

[WML20b] Wang, Jie, Victor Magron, and Jean-Bernard Lasserre.
*Chordal-TSSOS: a moment-SOS hierarchy that exploits term sparsity with chordal extension*.
arXiv preprint arXiv:2003.03210 (2020).
"""
struct Monomial{C<:CEG.AbstractCompletion} <: Pattern
    completion::C
    k::Int
    use_all_monomials::Bool
    function Monomial(
        completion::CEG.AbstractCompletion=CEG.ClusterCompletion(),
        k::Int=0,
        use_all_monomials::Bool=false,
    )
        return new{typeof(completion)}(completion, k, use_all_monomials)
    end
end

# Note: For some implementation of MultivariatePolynomials such as
# https://github.com/blegat/CondensedMatterSOS.jl,
# the product of monomials may be a term so wrap any multiplication of monomials by
# `MP.monomial`.

const MP = SumOfSquares.MP
function monomial_sparsity_graph(monos, P, use_all_monomials::Bool)
    g = CEG.LabelledGraph{eltype(monos)}()
    if use_all_monomials
        squares = nothing # No need to record this
    else
        squares = Set{eltype(monos)}()
    end
    for mono in monos
        CEG.add_node!(g, mono)
    end
    for a in monos
        for b in monos
            if MP.monomial(a * b) in P
                if a != b
                    CEG.add_edge!(g, a, b)
                elseif squares !== nothing
                    push!(squares, a)
                end
            end
        end
    end
    return g, squares
end
function _add_monos(add, monos, H, squares)
    for a in monos
        neighbors = H.graph.neighbors[H.n2int[a]]
        # If `!isempty!(neighbors)`, the monomial part of a clique with another
        # monomial so the diagonal entry should be nonzero.
        # If `a in squares`, the monomial square is needed so it is added.
        # e.g. if the polynomial is x^2 but the basis is `[x, 1]`,
        # `squares` is `Set([x])` so `x^2` will be added but not `1^2`.
        if squares === nothing || !isempty(neighbors) || a in squares
            add(a^2)
        end
        for bi in neighbors
            b = H.int2n[bi]
            add(MP.monomial(a * b))
        end
    end
end
function completion_with_squares(g, squares, completion)
    H, cliques = CEG.completion(g, completion)
    if squares !== nothing
        cliques = filter(monos -> length(monos) > 1 || (monos[1] in squares), collect(cliques))
    end
    return H, cliques
end
function monomial_sparsity_iteration(P, completion, use_all_monomials::Bool, monos)
    P_next = Set{eltype(P)}()
    g, squares = monomial_sparsity_graph(monos, P, use_all_monomials)
    H, cliques = completion_with_squares(g, squares, completion)
    _add_monos(mono -> push!(P_next, mono), monos, H, squares)
    return P_next, cliques
end
struct GeneratorP{PT, GT}
    P::PT
    generator_monos::GT
end
function Base.in(mono, g::GeneratorP)
    return any(g.generator_monos) do g_mono
        MP.monomial(mono * g_mono) in g.P
    end
end
function monomial_sparsity_iteration(P, completion, use_all_monomials::Bool, monos, multiplier_generator_monos)
    P_next, cliques = monomial_sparsity_iteration(P, completion, use_all_monomials, monos)
    multiplier_cliques = map(multiplier_generator_monos) do m
        multiplier_monos, generator_monos = m
        g, squares = monomial_sparsity_graph(multiplier_monos, GeneratorP(P, generator_monos), use_all_monomials)
        H, _cliques = completion_with_squares(g, squares, completion)
        _add_monos(multiplier_monos, H, squares) do mono
            for b in generator_monos
                push!(P_next, MP.monomial(mono * b))
            end
        end
        return _cliques
    end
    return P_next, (cliques, multiplier_cliques)
end
_monovec(cliques::AbstractVector{<:MP.AbstractMonomial}) = MP.monovec(cliques)
_monovec(cliques) = _monovec.(cliques)
function sparsity(monos::AbstractVector{<:MP.AbstractMonomial}, sp::Monomial,
                  gram_monos::AbstractVector = SumOfSquares.Certificate.monomials_half_newton_polytope(monos, tuple()),
                  args...)
    P = Set(monos)
    if sp.use_all_monomials
        for mono in gram_monos
            push!(P, mono^2)
        end
    end
    cliques = nothing
    iter = 0
    while iter < sp.k || iszero(sp.k)
        P_prev = P
        P, cliques = monomial_sparsity_iteration(P_prev, sp.completion, sp.use_all_monomials, gram_monos, args...)
        if iszero(iter)
            # If gram_monos + gram_monos !⊆ monos, then it's possible that P_prev !⊆ P
            P == P_prev && break
        else
            length(P) >= length(P_prev) || error("Set of monomials should be increasing in monomial sparsity iterations.")
            length(P) == length(P_prev) && break
        end
        iter += 1
    end
    return _monovec(cliques)
end
# This also checks that it is indeed a monomial basis
_monos(basis::MB.MonomialBasis) = basis.monomials
function _gram_monos(vars, certificate::SumOfSquares.Certificate.MaxDegree{CT, MB.MonomialBasis}) where CT
    return _monos(SumOfSquares.Certificate.maxdegree_gram_basis(MB.MonomialBasis, vars, certificate.maxdegree))
end
# poly = s0 + sum si gi
# where `s1` are the multipliers with basis `multiplier_generator_monos`
# we want to get the gram basis of `s0`
function _ideal_monos(poly_monos, multiplier_gram_monos)
    monos_set = Set(poly_monos)
    for (gram_monos, gen_monos) in multiplier_gram_monos
        for a in gram_monos
            for b in gram_monos
                for m in gen_monos
                    push!(monos_set, a * b * m)
                end
            end
        end
    end
    return MP.monovec(collect(monos_set))
end
# The ideal certificate should only ask for `MP.monomial`
struct DummyPolynomial{M}
    monomials::M
end
MP.monomials(p::DummyPolynomial) = p.monomials
MP.variables(p::DummyPolynomial) = MP.variables(p.monomials)
function sparsity(poly::MP.AbstractPolynomial, domain::SemialgebraicSets.BasicSemialgebraicSet, sp::Monomial, certificate::SumOfSquares.Certificate.AbstractPreorderCertificate)
    processed = SumOfSquares.Certificate.preprocessed_domain(certificate, domain, poly)
    multiplier_generator_monos = [
        (_monos(SumOfSquares.Certificate.multiplier_basis(certificate, index, processed)),
         MP.monomials(SumOfSquares.Certificate.generator(certificate, index, processed)))
        for index in SumOfSquares.Certificate.preorder_indices(certificate, processed)
    ]
    gram_monos = _monos(SumOfSquares.Certificate.gram_basis(
        SumOfSquares.Certificate.ideal_certificate(certificate),
        DummyPolynomial(_ideal_monos(MP.monomials(poly), multiplier_generator_monos)),
    ))
    cliques, multiplier_cliques = sparsity(MP.monomials(poly), sp, gram_monos, multiplier_generator_monos)
    return MB.MonomialBasis.(cliques), [MB.MonomialBasis.(clique) for clique in multiplier_cliques]
end
