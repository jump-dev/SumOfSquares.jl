struct MonomialSparsity <: Sparsity
    k::Int
end

const MP = SumOfSquares.MP
const CEG = SumOfSquares.Certificate.ChordalExtensionGraph
function monomial_sparsity_graph(monos, P)
    g = CEG.LabelledGraph{eltype(monos)}()
    for mono in monos
        CEG.add_node!(g, mono)
    end
    for a in monos
        for b in monos
            if a != b && (a * b) in P
                CEG.add_edge!(g, a, b)
            end
        end
    end
    return g
end
function monomial_sparsity_iteration(monos, P)
    g = monomial_sparsity_graph(monos, P)
    H, cliques = CEG.chordal_extension(g, CEG.GreedyFillIn())
    P_next = Set{eltype(P)}()
    for a in monos
        push!(P_next, a^2)
        for bi in H.graph.neighbors[H.n2int[a]]
            b = H.int2n[bi]
            push!(P_next, a * b)
        end
    end
    return P_next, cliques
end
function sparsity(monos::AbstractVector{<:MP.AbstractMonomial}, sp::MonomialSparsity)
    half_monos = SumOfSquares.Certificate.monomials_half_newton_polytope(monos, tuple())
    P = Set(monos)
    for mono in half_monos
        push!(P, mono^2)
    end
    cliques = nothing
    iter = 0
    while iter < sp.k || iszero(iter)
        P_prev = P
        P, cliques = monomial_sparsity_iteration(half_monos, P_prev)
        length(P) >= length(P_prev) || error("Set of monomials should be increasing in monomial sparsity iterations.")
        length(P) == length(P_prev) && break
    end
    return cliques
end
