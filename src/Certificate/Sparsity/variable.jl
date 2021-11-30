"""
    struct Variable <: Sparsity end

Variable or correlative sparsity as developed in [WSMM06].

[WSMM06] Waki, Hayato, Sunyoung Kim, Masakazu Kojima, and Masakazu Muramatsu. "Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity." SIAM Journal on Optimization 17, no. 1 (2006): 218-242.
"""
struct Variable <: Sparsity end

const CEG = ChordalExtensionGraph

function csp_graph(poly::MP.APL, ::FullSpace)
    G = CEG.LabelledGraph{MP.variable_union_type(poly)}()
    for mono in MP.monomials(poly)
        CEG.add_clique!(G, MP.effective_variables(mono))
    end
    return G
end

function csp_graph(poly::MP.APL, domain::AbstractAlgebraicSet)
    G = csp_graph(poly, FullSpace())
    for p in equalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function csp_graph(poly::MP.APL, domain::BasicSemialgebraicSet)
    G = csp_graph(poly, domain.V)
    for p in inequalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function chordal_csp_graph(poly::MP.APL, domain::AbstractBasicSemialgebraicSet)
    H, cliques = CEG.chordal_extension(csp_graph(poly, domain), CEG.GreedyFillIn())
    for clique in cliques
        sort!(clique, rev=true)
        unique!(clique)
    end
    return H, cliques

end

function sparsity(poly::MP.AbstractPolynomial, domain::BasicSemialgebraicSet,
                  sp::Variable, certificate::SumOfSquares.Certificate.Putinar)
    H, cliques = chordal_csp_graph(poly, domain)
    function bases(q)
        return [
            maxdegree_gram_basis(certificate.basis, clique,
                                 multiplier_maxdegree(certificate.maxdegree, q))
            for clique in cliques if variables(q) ⊆ clique
        ]
    end
    return bases(poly), map(bases, domain.p)
end
