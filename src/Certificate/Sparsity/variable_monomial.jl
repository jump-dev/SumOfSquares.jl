"""
    Sparsity.VariableMonomial(term = Sparsity.Monomial(ChordalCompletion(), 1, true))

Combined correlative and term sparsity (CS-TSSOS). Compute a minimum-fill
chordal extension of the variable graph, then apply `term` to the Gram and
localizing bases within those cliques. Each term-sparsity iteration uses
the support accumulated across all cliques. The default takes one iteration.

Equality localizers contribute to the support. Coefficient matching uses
the same ideal as the other sparsity patterns, so a truncated domain retains
all equality multiples allowed by its truncation degree.

See Wang, Magron and Lasserre, [Certifying global optimality of AC-OPF
solutions via sparse polynomial optimization](https://arxiv.org/abs/2109.10005),
Sections 3.1 and 3.2.
"""
struct VariableMonomial{T<:Monomial} <: Pattern
    term::T
end
VariableMonomial() = VariableMonomial(Monomial(ChordalCompletion(), 1, true))

# Within each degree, start with the earliest variables. Vertex ordering
# breaks ties during minimum fill, even at relaxation order one.
function _clique_monomials(clique, degree)
    return [
        m for d in 0:degree for m in Iterators.reverse(MP.monomials(clique, d))
    ]
end

# The graph for each localizer has its own vertices but shares its support
# with all the other graphs. Store a scalar localizer only once even if its
# polynomial spans several cliques.
function _localizers!(graphs, cliques, q, degree)
    start = length(graphs) + 1
    vars = MP.effective_variables(q)
    if degree >= 0 && !iszero(q)
        hosts = degree < 2 ? (vars,) : filter(c -> vars ⊆ c, cliques)
        for clique in hosts
            push!(
                graphs,
                (_clique_monomials(clique, fld(degree, 2)), MP.monomials(q)),
            )
        end
    end
    return start:length(graphs)
end

function _combined_bases(poly, domain, sp, maxdegree)
    basis = MB.explicit_basis(poly)
    vars = MP.variables(basis)
    if !(domain isa FullSpace)
        vars = sort!(union(vars, MP.variables(domain)); rev = true)
    end
    if isempty(vars)
        cliques = [vars]
    else
        _, cliques = chordal_csp_graph(
            basis,
            domain,
            maxdegree;
            variables = vars,
            algorithm = CEG.CliqueTrees.MinimalChordal(CEG.GreedyFillIn()),
        )
    end
    unit = MP.monomials(vars, 0:0)
    graphs = [(_clique_monomials(c, fld(maxdegree, 2)), unit) for c in cliques]
    support = Set(MB.keys_as_monomials(basis))
    gs = domain isa BasicSemialgebraicSet ? inequalities(domain) : ()
    multipliers = map(gs) do q
        union!(support, MP.monomials(q))
        return _localizers!(graphs, cliques, q, maxdegree - MP.maxdegree(q))
    end
    V = domain isa BasicSemialgebraicSet ? domain.V : domain
    for h in _equalities(V)
        union!(support, MP.monomials(h))
        if V isa SemialgebraicSets.TruncatedAlgebraicSet
            _localizers!(graphs, cliques, h, V.maxdegree - MP.maxdegree(h))
        end
    end
    blocks = _combined_term_blocks(graphs, sp.term, support)
    function bases(indices)
        result = typeof(basis)[]
        for i in indices, block in blocks[i]
            sub = MB.SubBasis{MB.Monomial}(MP.monomial_vector(block))
            ambient, _ = SA.promote_bases(sub, MB.implicit_basis(basis))
            push!(result, ambient)
        end
        return result
    end
    return bases(eachindex(cliques)), [bases(r) for r in multipliers]
end

function _combined_term_blocks(graphs, term, support)
    if term.use_all_monomials
        for (monos, weight) in graphs, m in monos, w in weight
            push!(support, MP.monomial(m^2 * w))
        end
    end
    step = 0
    while true
        step += 1
        next = empty(support)
        blocks = map(graphs) do (monos, weight)
            return _localizing_blocks!(
                next,
                support,
                term.completion,
                term.use_all_monomials,
                monos,
                weight,
            )
        end
        if step == term.k || next == support
            return blocks
        end
        support = next
    end
end

function sparsity(
    poly,
    sp::VariableMonomial,
    cert::SumOfSquares.Certificate.MaxDegree,
    domain = FullSpace(),
)
    return first(_combined_bases(poly, domain, sp, cert.maxdegree))
end

function sparsity(
    poly,
    domain::BasicSemialgebraicSet,
    sp::VariableMonomial,
    cert::SumOfSquares.Certificate.Putinar,
)
    return _combined_bases(poly, domain, sp, cert.maxdegree)
end
