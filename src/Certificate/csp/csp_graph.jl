export csp_graph, effective_variables, chordal_csp_graph

function effective_variables(p::APL)
    vs = [v for v in MP.variables(p)]
    ev = typeof(vs)()
    for v in vs
        if MP.maxdegree(p, v) > 0
            push!(ev, v)
        end
    end
    return ev
end

function csp_ojective_cliques(p::APL)
    return map(effective_variables, MP.monomials(p))
end

function csp_constraint_clique(p::APL)
    return effective_variables(p)
end

function csp_graph(::Type{T}, objective::APL, constraints::Vector{<:APL}, G::CEG.LabelledGraph{T} = CEG.LabelledGraph{T}() ) where T
    CEG.add_clique!.(G, csp_ojective_cliques(objective))
    CEG.add_clique!.(G, csp_constraint_clique.(constraints))
    return G
end

function csp_graph(objective::APL, constraints::Vector{<:APL}, G::CEG.LabelledGraph = CEG.LabelledGraph{MP.variable_union_type(objective)}() )
    CEG.add_clique!.(G, csp_ojective_cliques(objective))
    CEG.add_clique!.(G, csp_constraint_clique.(constraints))
    return G
end

function chordal_csp_graph(::Type{T}, objective::APL, constraints::Vector{<:APL}, G::CEG.LabelledGraph{T} = CEG.LabelledGraph{T}() ) where T
    return CEG.chordal_extension(csp_graph(T, objective, constraints, G), CEG.GreedyFillIn())
end

function chordal_csp_graph(objective::APL, constraints::Vector{<:APL}, G::CEG.LabelledGraph = CEG.LabelledGraph{MP.variable_union_type(objective)}() )
    return CEG.chordal_extension(csp_graph(objective, constraints, G), CEG.GreedyFillIn())
end

function chordal_csp_graph(objective::APL, G::CEG.LabelledGraph = CEG.LabelledGraph{MP.variable_union_type(objective)}() )
    return CEG.chordal_extension(csp_graph(objective, typeof(objective)[], G), CEG.GreedyFillIn())
end

function chordal_csp_graph(objective::APL, K::AbstractBasicSemialgebraicSet, G::CEG.LabelledGraph = CEG.LabelledGraph{MP.variable_union_type(objective)}() )
    constraints = Vector{typeof(objective)}()
    if K isa BasicSemialgebraicSet
        append!(constraints, inequalities(K))
    end
    if !(K isa FullSpace)
        append!(constraints, equalities(K))
    end
    return chordal_csp_graph(objective, constraints, G)
end
