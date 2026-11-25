using DynamicPolynomials
using PermutationGroups
using SumOfSquares
import Clarabel

function subgraph_counts(d)
    n = length(d)
    edges = (1 // 2) * sum(d)
    paths = (1 // 2) * sum(k -> k * (k - 1), d)
    single = (n - 2) * edges - 2 * paths
    empty = binomial(n, 3) - single - paths
    return edges, [empty, single, paths]
end

function mantel_bound(n)
    @polyvar d[1:n]
    edges, triples = subgraph_counts(d)
    group = PermGroup(perm"(1,2)", Perm([2:n; 1]))
    pattern = Symmetry.Pattern(group, Symmetry.VariablePermutation())
    model = Model(Clarabel.Optimizer)
    set_silent(model)
    @variable(model, β)
    @variable(model, λ[1:3] >= 0)
    @objective(model, Min, β)
    cert = @constraint(
        model,
        β - edges - sum(λ .* triples) in SOSCone(),
        symmetry = pattern,
    )
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return (; bound = value(β), weights = value.(λ), gram = gram_matrix(cert))
end

five = mantel_bound(5)
five.gram

solutions = [five, mantel_bound(6)]
map(5:6, solutions) do n, sol
    orders = sort([size(block.Q, 1) for block in sol.gram.blocks])
    return (; n, sol.bound, sol.weights, block_orders = orders)
end

n = 5
@polyvar d[1:n]
edges, triples = subgraph_counts(d)
squares = sum(k -> (k - n // 2)^2, d)
verified = n * (n^2 // 4 - edges) == triples[2] + squares

flag_difference = [-1, 1]
flag_difference * flag_difference' / 2

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
