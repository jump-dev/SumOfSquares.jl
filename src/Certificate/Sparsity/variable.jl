"""
    struct Sparsity.Variable <: Sparsity.Pattern end

Variable or correlative sparsity as developed in [WSMM06].

[WSMM06] Waki, Hayato, Sunyoung Kim, Masakazu Kojima, and Masakazu Muramatsu. "Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity." SIAM Journal on Optimization 17, no. 1 (2006): 218-242.
"""
struct Variable <: Pattern end

const CEG = ChordalExtensionGraph

function csp_graph(basis::MB.SubBasis{MB.Monomial}, ::FullSpace)
    G = CEG.LabelledGraph{MP.variable_union_type(eltype(basis.monomials))}()
    for mono in basis.monomials
        CEG.add_clique!(G, MP.effective_variables(mono))
    end
    return G
end

function csp_graph(basis::MB.SubBasis, domain::AbstractAlgebraicSet)
    G = csp_graph(basis, FullSpace())
    for p in equalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function csp_graph(basis::MB.SubBasis, domain::BasicSemialgebraicSet)
    G = csp_graph(poly, domain.V)
    for p in inequalities(domain)
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return G
end

function chordal_csp_graph(
    basis::MB.SubBasis,
    domain::AbstractBasicSemialgebraicSet,
)
    H, cliques =
        CEG.chordal_extension(csp_graph(basis, domain), CEG.GreedyFillIn())
    for clique in cliques
        sort!(clique, rev = true)
        unique!(clique)
    end
    return H, cliques
end

function sparsity(
    basis::MB.SubBasis{MB.Monomial},
    domain::BasicSemialgebraicSet,
    ::Variable,
    certificate::SumOfSquares.Certificate.Putinar,
)
    H, cliques = chordal_csp_graph(basis, domain)
    function bases(q)
        return [
            SumOfSquares.Certificate.maxdegree_gram_basis(
                certificate.multipliers_certificate.basis,
                clique,
                SumOfSquares.Certificate.multiplier_maxdegree(
                    certificate.maxdegree,
                    q,
                ),
            ) for clique in cliques if MP.variables(q) ⊆ clique
        ]
    end
    return bases(basis), map(bases, domain.p)
end
