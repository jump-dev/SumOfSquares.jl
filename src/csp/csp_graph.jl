include("ChordalExtensionGraph.jl")
const CEG = ChordalExtensionGraph

using MultivariatePolynomials
const MP = MultivariatePolynomials
const APL = AbstractPolynomialLike

function MP.degree(p::APL, v::AbstractVariable)
    d = 0
    for t in terms(p)
        d = max(d, MP.degree(t, v))
    end
    return d
end

function effective_variables(p::APL)
    vs = [v for v in variables(p)]
    ev = typeof(vs)()
    for v in vs
        if MP.degree(p, v) > 0
            push!(ev, v)
        end
    end
    return ev
end

function csp_ojective_cliques(p::APL)
    mv = monomials(p)
    return [effective_variables(m) for m in mv]
end

function csp_constraint_clique(p::APL)
    return effective_variables(p)
end

function csp_graph(::Type{T}, objective::APL, constraints::Vector{<:APL}, G::CEG.Graph{T} = CEG.Graph{T}() ) where T
    CEG.add_clique!.(G, csp_ojective_cliques(objective))
    CEG.add_clique!.(G, csp_constraint_clique.(constraints))
    return G
end

function chordal_csp_graph(::Type{T}, objective::APL, constraints::Vector{<:APL}, G::CEG.Graph{T} = CEG.Graph{T}() ) where T
    return CEG.chordal_extension(csp_graph(T, objective, constraints, G))
end

element_type(v::Array{T}) where T = T

function variable_type(v::Vector{<:APL})
    vararr = [i for i in variables(sum(vi for vi in v))]
    return element_type(vararr)
end

function variable_type(p::APL)
    return variable_type([p])
end

function csp_graph(objective::APL, constraints::Vector{<:APL}, G::CEG.Graph = CEG.Graph{variable_type(objective)}() )
    CEG.add_clique!.(G, csp_ojective_cliques(objective))
    CEG.add_clique!.(G, csp_constraint_clique.(constraints))
    return G
end

function chordal_csp_graph(objective::APL, constraints::Vector{<:APL}, G::CEG.Graph = CEG.Graph{variable_type(objective)}() )
    return CEG.chordal_extension(csp_graph(objective, constraints, G))
end

