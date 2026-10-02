using Test
import Random
using SumOfSquares
import MultivariateBases as MB
using DynamicPolynomials

function _moments(vars, centers, weights; maxdegree = 2)
    η = AtomicMeasure(vars, WeightedDiracMeasure.(centers, weights))
    basis = MB.SubBasis{MB.Monomial}(monomials(vars, 0:maxdegree))
    return moment_vector(η, basis)
end

@testset "heuristic_projection" begin
    @polyvar x[1:2]
    f1 = 2x[1]^4 - 8x[1]^3 + 8x[1]^2 + 2
    f2 = 4x[1]^4 - 32x[1]^3 + 88x[1]^2 - 96x[1] + 36
    K = @set x[1] >= 0 &&
         x[1] <= 3 &&
         x[2] >= 0 &&
         x[2] <= 4 &&
         x[2] <= f1 &&
         x[2] <= f2
    # `x[1]` is projected onto `[0, 3]` then `x[2]` onto
    # `[0, min(4, f1(x[1]), f2(x[1]))]`.
    @test heuristic_projection(K, x, [3, 4]) == [3, 0]
    @test heuristic_projection(K, x, [-1, 10]) == [0, 2]
    @test heuristic_projection(K, x, [2, 1]) == [2, 1]
    @test heuristic_projection(K, x, [1, 10]) == [1, 0]
    # `f1(0.5) = 3.125` and `f2(0.5) = 6.25`
    @test heuristic_projection(K, x, [0.5, 10]) ≈ [0.5, 3.125]
    @polyvar a b
    # Both variables appear in every constraint so the projection needs to
    # fix one of them. Fixing `a = 1.5` first leads to an infeasible problem
    # in `b` so it then tries fixing `b = 0.2`.
    @test heuristic_projection(@set(a^2 + b^2 <= 1), [a, b], [1.5, 0.2]) ≈
          [√0.96, 0.2]
    @test heuristic_projection(@set(a^2 + b^2 == 1), [a, b], [1.5, 0.2]) ≈
          [√0.96, 0.2]
    @test heuristic_projection(@set(a^2 + b^2 <= 1), [a, b], [0.5, 0.2]) ==
          [0.5, 0.2]
    @test isnothing(
        heuristic_projection(
            @set(a^2 + b^2 <= 1),
            [a, b],
            [1.5, 0.2],
            max_guesses = 1,
        ),
    )
    # Fixing any of the two variables at its current value leads to an
    # infeasible problem so the heuristic gives up even if the set is not
    # empty.
    @test isnothing(
        heuristic_projection(
            @set(a * b == 1 && a + b <= 3),
            [a, b],
            [0.1, 0.2],
        ),
    )
    @test isnothing(heuristic_projection(@set(a^2 <= -1), [a, b], [1, 2]))
    # `b` does not appear in the constraints so it is not modified
    @test heuristic_projection(@set(a == 1), [a, b], [2, 3]) == [1, 3]
    @test heuristic_projection(FullSpace(), [a, b], [2, 3]) == [2, 3]
    @test heuristic_projection(@set(a^2 == 2), [a, b], [-1, 0]) ≈ [-√2, 0]
    # Double root
    @test heuristic_projection(@set((a - 1)^2 == 0), [a, b], [3, 0]) ≈ [1, 0]
    @test_throws ErrorException heuristic_projection(@set(a == 1), [b], [2])
    # Once `a` is fixed to zero, `a * b - 1` is constant
    @test isnothing(
        heuristic_projection(@set(a == 0 && a * b >= 1), [a, b], [1, 1]),
    )
end

@testset "round_solution" begin
    @polyvar x y
    vars = [x, y]
    # MAX-CUT on a graph with a single edge
    K = @set x^2 == 1 && y^2 == 1
    μ = _moments(vars, [[1, -1], [-1, 1]], [0.5, 0.5])
    @test only(rounding_candidates(μ, vars, FirstMomentRounding())) ≈ [0, 0]
    gaussian =
        GaussianRounding(num_samples = 10, rng = Random.MersenneTwister(0))
    candidates = rounding_candidates(μ, vars, gaussian)
    @test length(candidates) == 11
    @test candidates[1] ≈ [0, 0]
    # The covariance matrix is `[1 -1; -1 1]`
    for c in candidates
        @test c[1] ≈ -c[2]
    end
    # The first moment is projected to a point with `x * y == 1`
    @test round_solution(μ, K, x * y) in [[-1, -1], [1, 1]]
    # Each sample `(t, -t)` is projected to `(sign(t), -sign(t))`
    @test round_solution(μ, K, x * y, rounding = gaussian) in [[-1, 1], [1, -1]]
    # The moments can also be given as a moment matrix
    ν = moment_matrix(μ, monomials(vars, 0:1))
    @test round_solution(ν, K, x * y, rounding = gaussian) in [[-1, 1], [1, -1]]
    # Without objective, the first feasible point is returned
    @test round_solution(μ, K) in [[-1, -1], [1, 1]]
    @test isnothing(round_solution(μ, @set(x^2 <= -1)))
    # Weighted by the mass
    μ = _moments(vars, [[1, 2]], [2])
    @test round_solution(μ, FullSpace()) ≈ [1, 2]
    μ = _moments(vars, [[1, 2]], [0])
    @test_throws ErrorException round_solution(μ, FullSpace())
    μ = _moments(vars, [[1, 2]], [1], maxdegree = 1)
    @test round_solution(μ, FullSpace()) ≈ [1, 2]
    @test_throws ErrorException round_solution(
        μ,
        FullSpace(),
        rounding = gaussian,
    )
    μ = _moments(vars, [[1, 2]], [1], maxdegree = 0)
    @test_throws ErrorException round_solution(μ, FullSpace())
    μ = _moments(vars, [[1, 2]], [1])
    @test_throws ErrorException rounding_candidates(
        μ,
        [x],
        FirstMomentRounding(),
    )
end
