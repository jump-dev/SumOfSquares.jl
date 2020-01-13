module ChordalExtensionGraph

export LabelledGraph, add_node!, add_edge!, add_clique!, chordal_extension

# With a `Vector{Vector{Int}}` with unsorted neighbors, computing fill-in
# would be inefficient.
# With a `Vector{Vector{Int}}` with sorted neighbors, adding an edge would be
# inefficient.
# With `Vector{Set{Int}}` both adding edges and computing fill-in is efficient.

struct Graph
    neighbors::Vector{Set{Int}}
    disabled::BitSet
end
Graph() = Graph(Set{Int}[], BitSet())
Base.broadcastable(g::Graph) = Ref(g)
Base.copy(G::Graph) = Graph(deepcopy(G.neighbors), copy(G.disabled))
function add_node!(G::Graph)
    push!(G.neighbors, Set{Int}())
    return length(G.neighbors)
end
num_nodes(G::Graph) = length(G.neighbors)
function add_edge!(G::Graph, i::Int, j::Int)
    if !(j in G.neighbors[i])
        push!(G.neighbors[i], j)
        push!(G.neighbors[j], i)
    end
    return (i, j)

end
function add_clique!(G::Graph, nodes::Vector{Int})
    for i in 1:(length(nodes) - 1)
        for j in (i + 1):length(nodes)
            add_edge!(G, nodes[i], nodes[j])
        end
    end
end
disable_node!(G::Graph, node::Int) = push!(G.disabled, node)
is_enabled(G::Graph, node::Int) = !(node in G.disabled)
enable_all_nodes!(G::Graph) = empty!(G.disabled)

"""
    neighbors(G::Graph, node::Int}

Return neighbors of `node` in `G`.
"""
function neighbors(G::Graph, node::Int)
    return G.neighbors[node]
end

function _num_edges_subgraph(G::Graph, nodes::Union{Vector{Int}, Set{Int}}, node::Int)
    neighs = neighbors(G, node)
    return count(nodes) do node
        is_enabled(G, node) && node in neighs
    end
end
function num_edges_subgraph(G::Graph, nodes::Union{Vector{Int}, Set{Int}})
    return mapreduce(+, nodes; init = 0) do node
        is_enabled(G, node) ? _num_edges_subgraph(G, nodes, node) : 0
    end
end

function num_missing_edges_subgraph(G::Graph, nodes::Union{Vector{Int}, Set{Int}})
    n = count(node -> is_enabled(G, node), nodes)
    # A clique is a completely connected graph. As such it has n*(n-1)/2 undirected
    # or equivalently n*(n-1) directed edges.
    return div(n * (n - 1) - num_edges_subgraph(G, nodes), 2)
end

"""
    fill_in(G::Graph{T}, i::T}

Return number of edges that need to be added to make the neighbors of `i` a clique.
"""
function fill_in(G::Graph, node::Int)
    return num_missing_edges_subgraph(G, neighbors(G, node))
end

"""
    is_clique(G::Graph{T}, x::Vector{T})

Return a `Bool` indication whether `x` is a clique in `G`.
"""
function is_clique(G::Graph, nodes::Vector{Int})
    return iszero(num_missing_edges_subgraph(G, nodes))
end

struct FillInCache
    graph::Graph
    fill_in::Vector{Int}
end
FillInCache(graph::Graph) = FillInCache(graph, [fill_in(graph, i) for i in 1:num_nodes(graph)])
Base.copy(G::FillInCache) = FillInCache(copy(G.graph), copy(G.fill_in))

num_nodes(G::FillInCache) = num_nodes(G.graph)
neighbors(G::FillInCache, node::Int) = neighbors(G.graph, node)
function add_edge!(G::FillInCache, i::Int, j::Int)
    ni = neighbors(G, i)
    nj = neighbors(G, j)
    if i in nj
        @assert j in ni
        return
    end
    @assert !(j in ni)
    for node in ni
        @assert node != i
        @assert node != j
        if node in nj
            G.fill_in[node] -= 1
        end
    end
    G.fill_in[i] += count(k -> is_enabled(G.graph, k), ni) - _num_edges_subgraph(G.graph, ni, j)
    G.fill_in[j] += count(k -> is_enabled(G.graph, k), nj) - _num_edges_subgraph(G.graph, nj, i)
    add_edge!(G.graph, i, j)
end
fill_in(G::FillInCache, node::Int) = G.fill_in[node]
function disable_node!(G::FillInCache, node::Int)
    for neighbor in neighbors(G, node)
        nodes = neighbors(G, neighbor)
        G.fill_in[neighbor] -= (length(nodes) - 1) - _num_edges_subgraph(G.graph, nodes, node)
    end
    disable_node!(G.graph, node)
end
is_enabled(G::FillInCache, node::Int) = is_enabled(G.graph, node)

"""
    struct LabelledGraph{T}
        n2int::Dict{T,Int}
        int2n::Vector{T}
        edges::Vector{Vector{Int}}
    end

Type to represend a graph with nodes of type T.
"""
struct LabelledGraph{T}
    n2int::Dict{T, Int}
    int2n::Vector{T}
    graph::Graph
end

Base.broadcastable(g::LabelledGraph) = Ref(g)
Base.copy(G::LabelledGraph{T}) where T = LabelledGraph(copy(G.n2int), copy(G.int2n), copy(G.graph))

function LabelledGraph{T}() where T
    return LabelledGraph(Dict{T, Int}(), T[], Graph())
end

function LabelledGraph()
    error("`LabelledGraph()` constructor is not valid, you need to specify the type of nodes as type parameter, e.g. `LabelledGraph{Int}()`.")
end

function Base.show(io::IO, G::LabelledGraph{T}) where T
    println(io, "$T-LabelledGraph with")
    println(io, "Nodes:")
    for node in G.int2n
        print(io, "$node, ")
    end
    if !isempty(G.int2n)
        println(io, "")
    end
    println(io, "Edges:")
    for (x, i) in G.n2int
        for j in G.graph.neighbors[i]
            println(io, "( $x, $(G.int2n[j]) )")
        end
    end
end

"""
    add_node!(G::LabelledGraph{T}, i::T) where T

Add the node i to graph G.
If i is already a node of G, only return the reference.
"""
function add_node!(G::LabelledGraph{T}, i::T) where T
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
    add_edge!(G::LabelledGraph{T}, i::T, j::T) where T

Add the unweighted edge (i, j) to graph G. Duplicate edges are not taken into account.
"""
function add_edge!(G::LabelledGraph{T}, xi::T, xj::T) where T
    add_edge!(G.graph, add_node!(G, xi), add_node!(G, xj))
    return (xi, xj)
end

"""
    add_edge!(G::LabelledGraph{T}, e::Tuple{T,T}) where T

Add the unweighted edge e to graph G. Duplicate edges are not taken into account.
"""
function add_edge!(G::LabelledGraph{T}, e::Tuple{T,T}) where T
    add_edge!(G, first(e), last(e))
end

"""
    add_clique!(G::LabelledGraph{T}, x::Vector{T}) where T

Add all elements of x as nodes to G and add edges such that x is fully
connected in G.
"""
function add_clique!(G::LabelledGraph{T}, x::Vector{T}) where T
    add_clique!(G.graph, add_node!.(G, x))
end

abstract type AbstractGreedyAlgorithm end

struct GreedyFillIn <: AbstractGreedyAlgorithm end
cache(G::Graph, ::GreedyFillIn) = FillInCache(G)
heuristic_value(G::FillInCache, node::Int, ::GreedyFillIn) = fill_in(G, node)

function _greedy_triangulation!(G, algo::AbstractGreedyAlgorithm)
    # Computes elimination graph `G` using equation (6.1) of [VA15].
    # generate cliques that are potentially maximal in the resulting graph G
    candidate_cliques = Vector{Int}[]
    for i in 1:num_nodes(G)
        node = argmin(map(1:num_nodes(G)) do node
            if is_enabled(G, node)
                heuristic_value(G, node, algo)
            else
                typemax(Int)
            end
        end)

        # look at all its neighbors. In G we want the neighbors to form a clique.
        neighbor_nodes = collect(neighbors(G, node))

        # add neighbors and node as a potentially maximal clique
        candidate_clique = [neighbor for neighbor in neighbor_nodes if is_enabled(G, neighbor)]
        push!(candidate_clique, node)
        push!(candidate_cliques, candidate_clique)

        # add edges to G to make the new candidate_clique at least potentially a clique
        for i in eachindex(neighbor_nodes)
            if is_enabled(G, i)
                for j in (i + 1):length(neighbor_nodes)
                    if is_enabled(G, j)
                        add_edge!(G, neighbor_nodes[i], neighbor_nodes[j])
                    end
                end
            end
        end
        disable_node!(G, node)
    end
    return candidate_cliques
end

function chordal_extension(G::Graph, algo::AbstractGreedyAlgorithm)
    H = copy(G)

    candidate_cliques = _greedy_triangulation!(cache(H, algo), algo)
    enable_all_nodes!(H)

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

    if length(Set(Iterators.flatten(maximal_cliques))) != num_nodes(H)
        error("Maximal cliques do not cover all nodes.")
    end

    return H, maximal_cliques
end

"""
    chordal_extension(G::LabelledGraph{T})

Return a chordal extension of G and the corresponding maximal cliques.

The algoritm is Algorithm 3 in [BA10] with the GreedyFillIn heuristic of Table I.

[BA10] Bodlaender, Hans L., and Arie MCA Koster.
Treewidth computations I. Upper bounds.
Information and Computation 208, no. 3 (2010): 259-275.
Utrecht University, Utrecht, The Netherlands www.cs.uu.nl
"""
function chordal_extension(G::LabelledGraph{T}, algo::AbstractGreedyAlgorithm) where T
    H, cliques = chordal_extension(G.graph, algo)
    return LabelledGraph{T}(G.n2int, G.int2n, H), [[G.int2n[i] for i in clique] for clique in cliques]
end
end #module
