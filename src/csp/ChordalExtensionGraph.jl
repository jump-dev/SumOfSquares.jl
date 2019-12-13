module ChordalExtensionGraph

export Graph, add_node!, add_edge!, add_clique!, chordal_extension

struct _Graph
    edges::Vector{Vector{Int}}
end
_Graph() = _Graph(Vector{Int}[])
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

"""
    chordal_extension(G::Graph{T})

Return a chordal extension of G and the corresponding maximal cliques.

Algorithm 1 in Treewidth Computations I.Upper BoundsHans L. BodlaenderArie M. C. A. Koster
Technical Report UU-CS-2008-032 September 2008 Department of Information and Computing Sciences
Utrecht University, Utrecht, The Netherlands www.cs.uu.nl
"""
function chordal_extension(G::Graph{T}) where T
    H = copy(G)
    nodes = copy(H.int2n)

    # bring into perfect elimination order based on fill_in
    # (less fill-in, earlier elimination)
    fill_nodes = Dict(node => fill_in(H, node) for node in nodes)
    sort!(nodes, by = i -> fill_nodes[i])
    elimination_order = Dict(nodes[i] => i for i in 1:length(nodes))

    # generate cliques that are potentially maximal in the resulting graph H
    candidate_cliques = Vector{Vector{T}}()
    for node in nodes
        elimination_node = elimination_order[node]

        # look at all its neighbors. In H we want the neighbors to form a clique.
        neighbor_nodes = neighbors(H, node)

        # add neighbors and node as a potentially maximal clique
        candidate_clique = push!(copy(neighbor_nodes), node)
        # add edges to H to make the new candidate_clique at least potentially a clique
        while !isempty(neighbor_nodes)

            neighbor = popfirst!(neighbor_nodes)
            other_neighbors = copy(neighbor_nodes)
            non_neighbors = Vector{T}()
            while !isempty(other_neighbors)
                next_neighbor = popfirst!(other_neighbors)

                # add an edge between neighbor and next_neighbor if both have more fill-in
                if elimination_order[neighbor] > elimination_node
                    if elimination_order[next_neighbor] > elimination_node
                        add_edge!(H, neighbor, next_neighbor)
                    else
                        if !is_clique(H,[neighbor, next_neighbor])
                            push!(non_neighbors, next_neighbor)
                        end
                    end
                else
                    if !is_clique(H,[neighbor, next_neighbor])
                        push!(non_neighbors, neighbor)
                    end
                end
            end
            # remove neighbors that are missing an edge in candidate_clique
            setdiff!(candidate_clique, non_neighbors)
        end
        push!(candidate_cliques, candidate_clique)
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

    # tests
    control = collect(Iterators.flatten(maximal_cliques))
    unique!(sort!(control))
    if sort!(H.int2n) != control
        error("Maximal cliques do not cover all nodes.")
    end

    return H, maximal_cliques
end
end #module
