struct MonomialSparsity <: Sparsity
    k::Int
end

const MP = SumOfSquares.MP
const CEG = SumOfSquares.Certificate.ChordalExtensionGraph
function monomial_sparsity_graph(monos, P)
    g = CEG.LabelledGraph{eltype(monos)}()
    for mono in monos
        CEG.add_node!(g, mono)
    end
    for a in monos
        for b in monos
            if a != b && (a * b) in P
                CEG.add_edge!(g, a, b)
            end
        end
    end
    return g
end
function _add_monos(add, monos, H)
    for a in monos
        add(a^2)
        for bi in H.graph.neighbors[H.n2int[a]]
            b = H.int2n[bi]
            add(a * b)
        end
    end
end
function monomial_sparsity_iteration(P, monos)
    P_next = Set{eltype(P)}()
    g = monomial_sparsity_graph(monos, P)
    H, cliques = CEG.chordal_extension(g, CEG.GreedyFillIn())
    _add_monos(mono -> push!(P_next, mono), monos, H)
    return P_next, cliques
end
struct GeneratorP{PT, GT}
    P::PT
    generator_monos::GT
end
function Base.in(mono, g::GeneratorP)
    return any(g.generator_monos) do g_mono
        (mono * g_mono) in g.P
    end
end
function monomial_sparsity_iteration(P, monos, multiplier_generator_monos)
    P_next, cliques = monomial_sparsity_iteration(P, monos)
    multiplier_cliques = map(multiplier_generator_monos) do m
        multiplier_monos, generator_monos = m
        g = monomial_sparsity_graph(multiplier_monos, GeneratorP(P, generator_monos))
        H, _cliques = CEG.chordal_extension(g, CEG.GreedyFillIn())
        _add_monos(multiplier_monos, H) do mono
            for b in generator_monos
                push!(P_next, mono * b)
            end
        end
        return _cliques
    end
    return P_next, (cliques, multiplier_cliques)
end
function sparsity(monos::AbstractVector{<:MP.AbstractMonomial}, sp::MonomialSparsity,
                  gram_monos::AbstractVector = SumOfSquares.Certificate.monomials_half_newton_polytope(monos, tuple()),
                  args...)
    P = Set(monos)
    for mono in gram_monos
        push!(P, mono^2)
    end
    cliques = nothing
    iter = 0
    while iter < sp.k || iszero(iter)
        P_prev = P
        P, cliques = monomial_sparsity_iteration(P_prev, gram_monos, args...)
        if iszero(iter)
            # If gram_monos + gram_monos !⊆ monos, then it's possible that P_prev !⊆ P
            P == P_prev && break
        else
            length(P) >= length(P_prev) || error("Set of monomials should be increasing in monomial sparsity iterations.")
            length(P) == length(P_prev) && break
        end
    end
    return cliques
end
# This also checks that it is indeed a monomial basis
_monos(basis::MB.MonomialBasis) = basis.monomials
function _gram_monos(vars, certificate::MaxDegree{CT, MB.MonomialBasis}) where CT
    return _monos(maxdegree_gram_basis(MB.MonomialBasis, vars, certificate.maxdegree))
end
function sparsity(poly::MP.AbstractPolynomial, domain::BasicSemialgebraicSet, sp::MonomialSparsity, certificate::AbstractPreorderCertificate)
    gram_monos = _gram_monos(
        reduce((v, q) -> unique!(sort!([v..., variables(q)...], rev=true)),
                  domain.p, init = variables(poly)),
        get(certificate, IdealCertificate())
    )
    processed = get(certificate, PreprocessedDomain(), domain, poly)
    multiplier_generator_monos = [
        (_monos(get(certificate, MultiplierBasis(), index, processed)),
         monomials(get(certificate, Generator(), index, processed)))
        for index in get(certificate, PreorderIndices(), processed)
    ]
    cliques, multiplier_cliques = sparsity(monomials(poly), sp, gram_monos, multiplier_generator_monos)
    return MB.MonomialBasis.(cliques), [MB.MonomialBasis.(clique) for clique in multiplier_cliques]
end
