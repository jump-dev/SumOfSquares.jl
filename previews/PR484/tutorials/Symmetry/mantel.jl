# # Mantel's theorem by symmetry reduction

#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Symmetry/mantel.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Symmetry/mantel.ipynb)

# How many edges can a triangle-free graph on `n` vertices have? A complete
# bipartite graph with parts as equal as possible has `⌊n²/4⌋` edges. We will
# find a sum-of-squares certificate that this is the maximum. The SDP has just
# two Gram blocks, of orders 2 and 1, for every `n ≥ 3`.
#
# We use the two rooted flags from Example 2.7 and the symmetry reduction
# from Section 5 of
#md # [Raymond2018](@cite).
#!md # [Raymond, Singh, and Thomas (2018)](https://doi.org/10.5802/alco.5).
# Choose a root vertex `i` and draw a second vertex uniformly from the other
# `n - 1` vertices. The flag records whether the pair is an edge. The
# second vertex is a neighbor with probability `d[i] / (n - 1)` and a
# non-neighbor with the complementary probability. Here `d[i]` is the degree
# of the root, the sum of the edge indicators incident to `i`. Thus the span
# of these flags over all roots has basis `1, d[1], …, d[n]`. A dense Gram
# matrix in this basis has order `n + 1`.

using DynamicPolynomials
using PermutationGroups
using SumOfSquares
import Clarabel
using Test #src

# ## Count induced triples
#
# Write `Hₖ` for the number of vertex triples that induce exactly `k` edges.
# Triangle-freeness leaves just `H₀, H₁, H₂`. Every pair of neighbors of a
# vertex gives a two-edge path with that vertex as its unique center.
# Consequently `H₂ = ∑ᵢ d[i](d[i] - 1)/2`. Each edge belongs to `n - 2`
# triples, so counting an edge together with a third vertex gives
# `(n - 2)m = H₁ + 2H₂`, where `m` is the number of edges. These counts use
# distinct vertices and are exact for finite `n`.

function subgraph_counts(d)
    n = length(d)
    edges = (1 // 2) * sum(d)
    paths = (1 // 2) * sum(k -> k * (k - 1), d)
    single = (n - 2) * edges - 2 * paths
    empty = binomial(n, 3) - single - paths
    return edges, [empty, single, paths]
end

# ## Search for the certificate
#
# We minimize an upper bound `β` on `m`. The solver chooses three nonnegative
# weights so that `β - m - (λ₀H₀ + λ₁H₁ + λ₂H₂)` is SOS as a polynomial in
# the formal degree variables. At the degree sequence of any triangle-free
# graph, each `Hₖ` counts triples and is nonnegative, which gives `m ≤ β`.
#
# Relabeling the vertices permutes the degrees. A transposition and an
# `n`-cycle generate all such relabelings. We pass this action to SumOfSquares
# through `symmetry` so SymbolicWedderburn computes the Gram blocks.

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

# For five vertices we can inspect the matrices and their bases directly.

five = mantel_bound(5)
five.gram

# Why do the block orders stay fixed? The constant and the sum of the degrees
# are both invariant, so they share a 2-by-2 Gram block. The degree deviations
# from the mean span the standard representation of the symmetric group.
# It occurs once, so one scalar Gram block suffices for all `n - 1` directions.
# The two blocks have four free entries. Together with the three weights and
# `β`, this gives eight scalar variables. Coefficient matching uses one equation
# for each of the four monomial orbits
# represented by `1`, `d[i]`, `d[i]^2`, and `d[i] * d[j]` with `i ≠ j`.
# Building the polynomials and computing the decomposition still depend on `n`.

solutions = [five, mantel_bound(6)]
map(5:6, solutions) do n, sol
    orders = sort([size(block.Q, 1) for block in sol.gram.blocks])
    @test orders == [1, 2] #src
    @test sol.bound ≈ n^2 / 4 atol = 1e-6 #src
    @test sol.weights ≈ [0, 1 / n, 0] atol = 1e-4 #src
    return (; n, sol.bound, sol.weights, block_orders = orders)
end

# ## Read off an exact proof
#
# The numerical solutions suggest `β = n²/4` and a weight of `1/n` on `H₁`.
# Substituting the counts above and completing the square gives the identity
#
# ```math
# n\left(\frac{n^2}{4} - m\right)
#   = H_1 + \sum_{i=1}^n \left(d_i - \frac{n}{2}\right)^2.
# ```
#
# We can verify the polynomial identity in exact arithmetic, for example at
# `n = 5`.

n = 5
@polyvar d[1:n]
edges, triples = subgraph_counts(d)
squares = sum(k -> (k - n // 2)^2, d)
verified = n * (n^2 // 4 - edges) == triples[2] + squares
@test verified #src

# Both terms on the right are nonnegative on every triangle-free graph.
# Since the edge count is an integer, `m ≤ ⌊n²/4⌋`, attained by the complete
# bipartite graph above. Dividing by `binomial(n, 2)` gives an edge density at
# most `n / (2(n - 1))`, with sharp limiting density `1/2`.
#
# How does this recover the flag certificate in the paper? In density units,
# the coefficient of the one-edge triple density is
# `binomial(n, 3) / (n * binomial(n, 2)) = (n - 2)/(3n)`, which tends to
# `1/3`. Write `F₁ = d[i]/(n - 1)` and `F₀ = 1 - F₁` for the two flag
# probabilities at a uniformly random root. The square term in density units
# is the root average of `2(n - 1)/n * (F₁ - n/(2(n - 1)))²`. In the limit
# this gives the average of `(F₁ - F₀)²/2`. Expanding the square in the
# ordered flag basis `(F₀, F₁)` recovers the paper's Gram matrix.

flag_difference = [-1, 1]
flag_difference * flag_difference' / 2
