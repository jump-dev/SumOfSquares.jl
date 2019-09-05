__precompile__()

module ChordalExtensionGraph

using SparseArrays

export Graph, add_node!, add_edge!, add_clique!, sub_graph, chordal_extension, adjacency_matrix


"""
    struct Graph{T}
        nodes::Vector{T}
        edges::Vector{Tuple{T, T}}
    end

Type to represend a graph with nodes of type T.
"""
struct Graph{T}
    nodes::Vector{T}
    edges::Vector{Tuple{T, T}}
end

Base.broadcastable(g::Graph) = Ref(g)
Base.copy(G::Graph{T}) where T = Graph(copy(G.nodes), copy(G.edges))

function Graph{T}() where T
    return Graph(T[], Tuple{T, T}[])
end

function Graph()
    show("Specify type of nodes!")
end

function Base.show(io::IO, G::Graph{T}) where T
    println(io, "$T-Graph with")
    println(io, "Nodes:")
    for node in G.nodes
        print(io, "$node, ")
    end
    if !isempty(G.nodes)
        println("")
    end
    println(io, "Edges:")
    for edge in  G.edges
        println(io, edge)
    end
end

"""
    add_node!(G::Graph{T}, i::T) where T

Add the node i to graph G. 
If i is already a node of G, only return the reference.
"""
function add_node!(G::Graph{T}, i::T) where T
    if !(i in G.nodes)
        push!(G.nodes, i)
    end
    return i
end

"""
    add_edge!(G::Graph{T}, i::T, j::T) where T

Add the unweighted edge (i, j) to graph G. Duplicate edges are not taken into account. 
"""
function add_edge!(G::Graph{T}, xi::T, xj::T) where T
    i, j = add_node!.(G, (xi, xj))
    append!(G.edges, [(i, j), (j, i)])
    unique!(sort!(G.edges))
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
    n_edges = length(G.edges)
    while !isempty(x)
        y = pop!(x)
        for z in x
            add_edge!(G, y, z)
        end
    end
end

"""
    sub_graph(G::Graph{T}, x::Vector{T}) 

Return restriction of G to x.
"""
function sub_graph(G::Graph{T}, x::Vector{T}) where T
    edges = Tuple{T, T}[]
    for edge in G.edges
        if first(edge) in x && last(edge) in x
            push!(edges, edge)
        end
    end
    return Graph(x, sort!(edges))
end

"""
    neighbors(G::Graph{T}, i::T}

Return neighbors of i in G.
"""
function neighbors(G::Graph{T}, i::T) where T
    neighbor_nodes = Vector{T}()
    bool = true
    edges = copy(G.edges)
    while bool&&!isempty(edges)
        edge = popfirst!(edges)
        if  edge[1] == i
            push!(neighbor_nodes, edge[2])
        elseif  edge[1] > i
            bool = false
        end
    end
    return neighbor_nodes
end

"""
    fill_in(G::Graph{T}, i::T}

Return number of edges that need to be added to make the neighors of i a clique.
Also return the set of neighbors
"""
function fill_in(G::Graph{T}, i::T) where T
    N = neighbors(G, i)
    S = sub_graph(G, N)
    n = length(S.nodes)
    return Int((n)*(n-1)/2-length(S.edges)/2)
end

"""
    is_clique(G::Graph{T}, x::Vector{T})

Return true if x is a clique in G.
"""
function is_clique(G::Graph{T}, x::Vector{T}) where T
    S = sub_graph(G, x)
    n = length(x)
    # A clique is a completely connected graph. As such it has n*(n-1)/2 undirected
    # or equivalently n*(n-1) directed edges.
    return length(S.edges)==n*(n-1)
end

function contains(bigvector::Vector{T}, Smallvector::Vector{T}) where T
    smallvector = copy(Smallvector)
    bool = !isempty(bigvector)
    while bool && !isempty(smallvector)
        foo = pop!(smallvector)
        bool = findfirst(x->x==foo, bigvector) isa Int
    end
   return bool
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
    nodes = copy(H.nodes)

    # bring into perfect elimination order based on fill_in
    # (more fill-in, earlier elimination)
    fill_nodes = Dict(node => fill_in(H, node) for node in nodes)
    sort!(nodes, lt = (i, j) -> fill_nodes[i] < fill_nodes[j])
    elimination_order = Dict(nodes[i] => i for i in 1:length(nodes))
     
    # generate cliques that are potentially maximal in the resulting graph H
    candidate_cliques = Vector{Vector{T}}()    
    while !isempty(nodes)

        # take the node with the smallest fill-in
        node = popfirst!(nodes)
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
                        push!(non_neighbors, next_neighbor) 
                    end
                else
                    push!(non_neighbors, neighbor) 
                end
            end
            # remove neighbors that are missing an edge in candidate_clique
            setdiff!(candidate_clique, non_neighbors)
        end
        push!(candidate_cliques, candidate_clique)
        sort!(H.edges)
        unique!(H.edges)
    end
    sort!.(candidate_cliques) 
    unique!(candidate_cliques) 
    # check whether candidate cliques are actually cliques
    candidate_cliques = candidate_cliques[[is_clique(H, clique) for clique in candidate_cliques]] 
    sort!(candidate_cliques, lt = (x,y)-> length(x)<length(y))
    maximal_cliques = [pop!(candidate_cliques)]
    while !isempty(candidate_cliques)
        clique = pop!(candidate_cliques)
        bool = true
        cliques = copy(maximal_cliques)
        while bool && !isempty(cliques)
            other_clique = pop!(cliques)
            bool = !contains(other_clique, clique)
        end
        if bool
            push!(maximal_cliques, clique)
        end        
    end
    return H, maximal_cliques
end


"""
    adjacency_matrix(G::Graph{T}) where T

Create the adjacency matrix of G.
"""
function adjacency_matrix(G::Graph{T}) where T
    n = length(G.nodes)
    id = Dict(G.nodes[i] => i for i=1:n)
    I = [(i, i) for i=1:n]
    for edge in G.edges
        push!(I, (id[edge[1]], id[edge[2]]))
    end
    unique!(I)
    return sparse(first.(I), last.(I), ones(Int, length(I)))
end
end #module
