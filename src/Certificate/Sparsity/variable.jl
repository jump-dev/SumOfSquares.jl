"""
    struct Sparsity.Variable <: Sparsity.Pattern end

Variable or correlative sparsity as developed in [WSMM06].
For a scalar localizing constraint, the variables of each monomial must
belong to a common clique. Different monomials can belong to different
cliques. This applies to inequalities and to equalities in a truncated
ideal at the top relaxation degree [WML, Section 3.1].

[WML] Jie Wang, Victor Magron, and Jean-Bernard Lasserre.
[Certifying global optimality of AC-OPF solutions via sparse polynomial optimization](https://arxiv.org/abs/2109.10005).

[WSMM06] Waki, Hayato, Sunyoung Kim, Masakazu Kojima, and Masakazu Muramatsu. "Sums of squares and semidefinite program relaxations for polynomial optimization problems with structured sparsity." SIAM Journal on Optimization 17, no. 1 (2006): 218-242.
"""
struct Variable <: Pattern end

const CEG = ChordalExtensionGraph

# At the highest constraint degree the WML localizing matrix is scalar.
# Each term then needs a clique, but terms need not share a clique.
function _constraint_edges!(G, p, maxdegree)
    if maxdegree !== nothing && 0 <= maxdegree - MP.maxdegree(p) < 2
        for mono in MP.monomials(p)
            CEG.add_clique!(G, MP.effective_variables(mono))
        end
    else
        CEG.add_clique!(G, MP.effective_variables(p))
    end
    return
end

function csp_graph(
    basis::MB.SubBasis{MB.Monomial},
    ::FullSpace,
    maxdegree = nothing,
)
    G = CEG.LabelledGraph{MP.variable_union_type(MP.monomial_type(basis))}()
    for mono in MB.keys_as_monomials(basis)
        CEG.add_clique!(G, MP.effective_variables(mono))
    end
    return G
end

function csp_graph(
    basis::MB.SubBasis,
    domain::AbstractAlgebraicSet,
    maxdegree = nothing,
)
    G = csp_graph(basis, FullSpace())
    degree =
        domain isa SemialgebraicSets.TruncatedAlgebraicSet ? domain.maxdegree :
        nothing
    for p in equalities(domain)
        _constraint_edges!(G, p, degree)
    end
    return G
end

function csp_graph(
    basis::MB.SubBasis,
    domain::BasicSemialgebraicSet,
    maxdegree = nothing,
)
    G = csp_graph(basis, domain.V, maxdegree)
    for p in inequalities(domain)
        _constraint_edges!(G, p, maxdegree)
    end
    return G
end

function chordal_csp_graph(
    basis::MB.SubBasis,
    domain::AbstractBasicSemialgebraicSet,
    maxdegree = nothing,
)
    H, cliques = CEG.chordal_extension(
        csp_graph(basis, domain, maxdegree),
        CEG.GreedyFillIn(),
    )
    for clique in cliques
        sort!(clique, rev = true)
        unique!(clique)
    end
    return H, cliques
end

function sparsity(
    poly,
    domain::BasicSemialgebraicSet,
    ::Variable,
    certificate::SumOfSquares.Certificate.Putinar,
)
    basis = MB.explicit_basis(poly)
    _, cliques = chordal_csp_graph(basis, domain, certificate.maxdegree)
    function clique_bases(q, vars)
        degree = SumOfSquares.Certificate.multiplier_maxdegree(
            certificate.maxdegree,
            q,
        )
        if 0 <= degree < 2
            return [
                SumOfSquares.Certificate.maxdegree_gram_basis(
                    certificate.multipliers_certificate.gram_basis,
                    vars,
                    0,
                ),
            ]
        end
        return [
            SumOfSquares.Certificate.maxdegree_gram_basis(
                certificate.multipliers_certificate.gram_basis,
                clique,
                degree,
            ) for clique in cliques if vars ⊆ clique
        ]
    end
    ideal_bases = [
        _maxdegree_basis(
            MB.implicit_basis(certificate.ideal_certificate),
            clique,
            certificate.maxdegree,
        ) for clique in cliques
    ]
    preorder_bases =
        map(q -> clique_bases(q, MP.effective_variables(q)), domain.p)
    return ideal_bases, preorder_bases
end
