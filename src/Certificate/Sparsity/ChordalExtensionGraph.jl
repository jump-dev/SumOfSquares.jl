module ChordalExtensionGraph

import CliqueTrees
import DataStructures
import SparseArrays

export AbstractCompletion, ChordalCompletion, ClusterCompletion
export LabelledGraph, add_node!, add_edge!, add_clique!, chordal_extension

abstract type AbstractGreedyAlgorithm end

abstract type AbstractCompletion end

struct ClusterCompletion <: AbstractCompletion end

struct ChordalCompletion{A<:CliqueTrees.EliminationAlgorithm} <:
       AbstractCompletion
    algo::A
end

# the default algorithm is the greedy minimum fill heuristic
const GreedyFillIn = CliqueTrees.MF

function ChordalCompletion()
    return ChordalCompletion(GreedyFillIn())
end

# With a `Vector{Vector{Int}}` with unsorted neighbors, computing fill-in
# would be inefficient.
# With a `Vector{Vector{Int}}` with sorted neighbors, adding an edge would be
# inefficient.
# With `Vector{Set{Int}}` both adding edges and computing fill-in is efficient.
struct Graph
    neighbors::Vector{Set{Int}}
end
Graph() = Graph(Set{Int}[])
Base.broadcastable(g::Graph) = Ref(g)
Base.copy(G::Graph) = Graph(deepcopy(G.neighbors))
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
    n = length(nodes)

    for i in 1:(n-1)
        for j in (i+1):n
            add_edge!(G, nodes[i], nodes[j])
        end
    end
end

"""
    neighbors(G::Graph, node::Int}

Return neighbors of `node` in `G`.
"""
function neighbors(G::Graph, node::Int)
    return G.neighbors[node]
end

"""
    is_clique(G::Graph{T}, x::Vector{T})

Return a `Bool` indication whether `x` is a clique in `G`.
"""
function is_clique(G::Graph, nodes::Vector{Int})
    n = length(nodes)

    for i in 1:(n-1)
        for j in (i+1):n
            if nodes[j] ∉ neighbors(G, nodes[i])
                return false
            end
        end
    end

    return true
end

"""
    struct LabelledGraph{T}
        n2int::Dict{T,Int}
        int2n::Vector{T}
        graph::Graph
    end

Type to represend a graph with nodes of type T.
"""
struct LabelledGraph{T}
    n2int::Dict{T,Int}
    int2n::Vector{T}
    graph::Graph
end

Base.broadcastable(g::LabelledGraph) = Ref(g)
function Base.copy(G::LabelledGraph{T}) where {T}
    return LabelledGraph(copy(G.n2int), copy(G.int2n), copy(G.graph))
end

function LabelledGraph{T}() where {T}
    return LabelledGraph(Dict{T,Int}(), T[], Graph())
end

function LabelledGraph()
    return error(
        "`LabelledGraph()` constructor is not valid, you need to specify the type of nodes as type parameter, e.g. `LabelledGraph{Int}()`.",
    )
end

function Base.show(io::IO, G::LabelledGraph{T}) where {T}
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
function add_node!(G::LabelledGraph{T}, i::T) where {T}
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
function add_edge!(G::LabelledGraph{T}, xi::T, xj::T) where {T}
    add_edge!(G.graph, add_node!(G, xi), add_node!(G, xj))
    return (xi, xj)
end

"""
    add_edge!(G::LabelledGraph{T}, e::Tuple{T,T}) where T

Add the unweighted edge e to graph G. Duplicate edges are not taken into account.
"""
function add_edge!(G::LabelledGraph{T}, e::Tuple{T,T}) where {T}
    return add_edge!(G, first(e), last(e))
end

"""
    add_clique!(G::LabelledGraph{T}, x::Vector{T}) where T

Add all elements of x as nodes to G and add edges such that x is fully
connected in G.
"""
function add_clique!(G::LabelledGraph{T}, x::Vector{T}) where {T}
    return add_clique!(G.graph, add_node!.(G, x))
end

"""
    completion(G::Graph, comp::ClusterCompletion)

Return a cluster completion of `G` and the corresponding maximal cliques.
"""
function completion(G::Graph, ::ClusterCompletion)
    H = copy(G)
    union_find = DataStructures.IntDisjointSets(num_nodes(G))
    for from in 1:num_nodes(G)
        for to in neighbors(G, from)
            DataStructures.union!(union_find, from, to)
        end
    end
    cliques = [Int[] for i in 1:DataStructures.num_groups(union_find)]
    ids = zeros(Int, num_nodes(G))
    k = 0
    for node in 1:num_nodes(G)
        root = DataStructures.find_root!(union_find, node)
        if iszero(ids[root])
            k += 1
            ids[root] = k
        end
        push!(cliques[ids[root]], node)
    end
    @assert k == length(cliques)
    for clique in cliques
        add_clique!(H, clique)
    end
    return H, cliques
end

"""
    completion(G::Graph, comp::ChordalCompletion)

Return a chordal extension of `G` and the corresponding maximal cliques.

The algoritm is Algorithm 3 in [BA10] with the GreedyFillIn heuristic of Table I.

[BA10] Bodlaender, Hans L., and Arie MCA Koster.
Treewidth computations I. Upper bounds.
Information and Computation 208, no. 3 (2010): 259-275.
Utrecht University, Utrecht, The Netherlands www.cs.uu.nl
"""
function completion(G::Graph, comp::ChordalCompletion)
    # construct a copy H of the graph G
    H = copy(G)

    # construct adjacency matrix of H
    matrix = SparseArrays.spzeros(Bool, num_nodes(H), num_nodes(H))

    for j in axes(matrix, 2)
        list = neighbors(H, j)
        SparseArrays.getcolptr(matrix)[j+1] =
            SparseArrays.getcolptr(matrix)[j] + length(list)
        append!(SparseArrays.rowvals(matrix), list)
        append!(SparseArrays.nonzeros(matrix), zeros(Bool, length(list)))
    end

    # sort row indices of adjacency matrix
    matrix = copy(transpose(matrix))

    # construct tree decomposition of H
    label, tree = CliqueTrees.cliquetree(matrix; alg = comp.algo)

    # construct vector of cliques and append fill edges to H
    cliques = Vector{Vector{Int}}(undef, length(tree))

    for (i, v) in enumerate(tree)
        clique = cliques[i] = label[v]
        add_clique!(H, clique)
    end

    return H, cliques
end

function completion(G::LabelledGraph{T}, comp::AbstractCompletion) where {T}
    H, cliques = completion(G.graph, comp)
    return LabelledGraph{T}(G.n2int, G.int2n, H),
    [[G.int2n[i] for i in clique] for clique in cliques]
end

chordal_extension(G, algo) = completion(G, ChordalCompletion(algo))

end #module
