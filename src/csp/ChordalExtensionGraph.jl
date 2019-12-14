module ChordalExtensionGraph

export Graph, add_node!, add_edge!, add_clique!, chordal_extension

struct _Graph
    edges::Vector{Vector{Int}}
end
_Graph() = _Graph(Vector{Int}[])
Base.broadcastable(g::_Graph) = Ref(g)
Base.copy(G::_Graph) = _Graph(deepcopy(G.edges))
function add_node!(G::_Graph)
    push!(G.edges, Int[])
    return length(G.edges)
end
function add_edge!(G::_Graph, i::Int, j::Int)
    if !(j in G.edges[i])
        push!(G.edges[i], j)
        push!(G.edges[j], i)
    end
    return (i, j)

end
function add_clique!(G::_Graph, nodes::Vector{Int})
    for i in 1:(length(nodes) - 1)
        for j in (i + 1):length(nodes)
            add_edge!(G, nodes[i], nodes[j])
        end
    end
end
function neighbors(G::_Graph, i::Int)
    return G.edges[i]
end
function num_edges_subgraph(G::_Graph, nodes::Vector{Int})
    return sum(nodes) do node
        length(nodes ∩ neighbors(G, node))
    end
end
function num_missing_edges_subgraph(G::_Graph, nodes::Vector{Int})
    n = length(nodes)
    # A clique is a completely connected graph. As such it has n*(n-1)/2 undirected
    # or equivalently n*(n-1) directed edges.
    return div(n * (n - 1) - num_edges_subgraph(G, nodes), 2)
end
function fill_in(G::_Graph, node::Int)
    return num_missing_edges_subgraph(G, neighbors(G, node))
end
function is_clique(G::_Graph, nodes::Vector{Int})
    return iszero(num_missing_edges_subgraph(G, nodes))
end


"""
    struct Graph{T}
        n2int::Dict{T,Int}
        int2n::Vector{T}
        edges::Vector{Vector{Int}}
    end

Type to represend a graph with nodes of type T.
"""
struct Graph{T}
    n2int::Dict{T, Int}
    int2n::Vector{T}
    graph::_Graph
end

Base.broadcastable(g::Graph) = Ref(g)
Base.copy(G::Graph{T}) where T = Graph(copy(G.n2int), copy(G.int2n), copy(G.graph))

function Graph{T}() where T
    return Graph(Dict{T, Int}(), T[], _Graph())
end

function Graph()
    error("`Graph()` constructor is not valid, you need to specify the type of nodes as type parameter, e.g. `Graph{Int}()`.")
end

function Base.show(io::IO, G::Graph{T}) where T
    println(io, "$T-Graph with")
    println(io, "Nodes:")
    for node in G.int2n
        print(io, "$node, ")
    end
    if !isempty(G.int2n)
        println(io, "")
    end
    println(io, "Edges:")
    for (x, i) in G.n2int
        for j in G.graph.edges[i]
            println(io, "( $x, $(G.int2n[j]) )")
        end
    end
end

"""
    add_node!(G::Graph{T}, i::T) where T

Add the node i to graph G.
If i is already a node of G, only return the reference.
"""
function add_node!(G::Graph{T}, i::T) where T
    if haskey(G.n2int, i)
        idx = G.n2int[i]
    else
        push!(G.int2n, i)
        idx = add_node!(G.graph)
        @assert idx == length(G.int2n)
        G.n2int[i] = idx
    end
    return idx
end

"""
    add_edge!(G::Graph{T}, i::T, j::T) where T

Add the unweighted edge (i, j) to graph G. Duplicate edges are not taken into account.
"""
function add_edge!(G::Graph{T}, xi::T, xj::T) where T
    add_edge!(G.graph, add_node!(G, xi), add_node!(G, xj))
    return (xi, xj)
end

"""
    add_edge!(G::Graph{T}, e::Tuple{T,T}) where T

Add the unweighted edge e to graph G. Duplicate edges are not taken into account.
"""
function add_edge!(G::Graph{T}, e::Tuple{T,T}) where T
    add_edge!(G, first(e), last(e))
end

"""
    add_clique!(G::Graph{T}, x::Vector{T}) where T

Add all elements of x as nodes to G and add edges such that x is fully
connected in G.
"""
function add_clique!(G::Graph{T}, x::Vector{T}) where T
    add_clique!(G.graph, add_node!.(G, x))
end

"""
    neighbors(G::Graph{T}, i::T}

Return neighbors of i in G.
"""
function neighbors(G::Graph{T}, i::T) where T
    return [G.int2n[j] for j in neighbors(G.graph, G.n2int[i])]
end

function n_edges(G::Graph{T}) where T
    return sum(length(neighbor) for neighbor in G.graph.edges)
end

"""
    fill_in(G::Graph{T}, i::T}

Return number of edges that need to be added to make the neighbors of `i` a clique.
"""
function fill_in(G::Graph{T}, i::T) where T
    return fill_in(G.graph, G.n2int[i])
end

"""
    is_clique(G::Graph{T}, x::Vector{T})

Return a `Bool` indication whether `x` is a clique in `G`.
"""
function is_clique(G::Graph{T}, x::Vector{T}) where T
    return is_clique(G.graph, [G.n2int[i] for i in x])
end

function chordal_extension(G::_Graph)
    H = copy(G)
    num_nodes = length(H.edges)

    # Bring into perfect elimination order based on `fill_in`.
    # (less fill-in, earlier elimination)
    σ = sortperm(fill_in.(H, 1:num_nodes))
    elimination_order = zeros(Int, num_nodes)
    for i in 1:num_nodes
        elimination_order[σ[i]] = i
    end

    # Computes elimination graph `H` using equation (6.1) of [VA15].
    # generate cliques that are potentially maximal in the resulting graph H
    candidate_cliques = Vector{Int}[]
    for node in σ
        elimination_node = elimination_order[node]

        # look at all its neighbors. In H we want the neighbors to form a clique.
        neighbor_nodes = neighbors(H, node)

        # add neighbors and node as a potentially maximal clique
        candidate_clique = [neighbor for neighbor in neighbor_nodes if elimination_order[neighbor] > elimination_node]
        push!(candidate_clique, node)
        push!(candidate_cliques, candidate_clique)

        # add edges to H to make the new candidate_clique at least potentially a clique
        for i in eachindex(neighbor_nodes)
            if elimination_order[neighbor_nodes[i]] > elimination_node
                for j in (i + 1):length(neighbor_nodes)
                    if elimination_order[neighbor_nodes[j]] > elimination_node
                        add_edge!(H, neighbor_nodes[i], neighbor_nodes[j])
                    end
                end
            end
        end
    end

    sort!.(candidate_cliques)
    unique!(candidate_cliques)

    # check whether candidate cliques are actually cliques
    candidate_cliques = candidate_cliques[[is_clique(H, clique) for clique in candidate_cliques]]
    sort!(candidate_cliques, by = x -> length(x))
    reverse!(candidate_cliques) # TODO use `rev=true` in `sort!`.
    maximal_cliques = [first(candidate_cliques)]
    for clique in Iterators.drop(candidate_cliques, 1)
        if all(other_clique -> !(clique ⊆ other_clique), maximal_cliques)
            push!(maximal_cliques, clique)
        end
    end

    if length(Set(Iterators.flatten(maximal_cliques))) != num_nodes
        error("Maximal cliques do not cover all nodes.")
    end

    return H, maximal_cliques
end

"""
    chordal_extension(G::Graph{T})

Return a chordal extension of G and the corresponding maximal cliques.

The algoritm is Algorithm 3 in [BA10] with the GreedyFillIn heuristic of Table I.

[BA10] Bodlaender, Hans L., and Arie MCA Koster.
Treewidth computations I. Upper bounds.
Information and Computation 208, no. 3 (2010): 259-275.
Utrecht University, Utrecht, The Netherlands www.cs.uu.nl
"""
function chordal_extension(G::Graph{T}) where T
    H, cliques = chordal_extension(G.graph)
    return Graph{T}(G.n2int, G.int2n, H), [[G.int2n[i] for i in clique] for clique in cliques]
end
end #module
