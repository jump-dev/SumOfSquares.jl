struct VariableSparsity <: Sparsity end

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

