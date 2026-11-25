module TestOptimizer

using Test

using LinearAlgebra
import MathOptInterface as MOI
using JuMP
using SumOfSquares
using DynamicPolynomials
import MultivariateBases as MB

import Clarabel
const SOLVER =
    optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

function test_polynomial_bases()
    @polyvar x
    for B in (MB.Monomial, MB.ScaledMonomial, MB.Chebyshev)
        p = MB.convert_basis(MB.FullBasis{B}(x), polynomial(1 + x^2))
        model = SOSModel(SOLVER)
        @variable(model, t)
        cref = @constraint(model, p - t in SOSCone())
        @objective(model, Max, t)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(t) ≈ 1 atol = 1e-6
        @test polynomial(gram_matrix(cref)) ≈ p - value(t) atol = 1e-6
        @test polynomial(sos_decomposition(cref)) ≈ p - value(t) atol = 1e-6
    end
    return
end

function test_noncommutative_polynomial()
    @ncpolyvar x y
    q = x * y + x^2
    p = q' * q
    model = SOSModel(SOLVER)
    cref = @constraint(model, p in SOSCone())
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test polynomial(gram_matrix(cref)) ≈ p atol = 1e-6
    @test polynomial(sos_decomposition(cref)) ≈ p atol = 1e-6
    return
end

function test_optimizer_attributes()
    optimizer = SumOfSquares.Optimizer(SOLVER)
    @test optimizer isa SumOfSquares.Optimizer{Float64}
    @test MOI.get(optimizer, MOI.SolverName()) == "SumOfSquares"
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
    @test MOI.get(optimizer, MOI.ResultCount()) == 0
    list = MOI.get(optimizer, MOI.Bridges.ListOfNonstandardBridges{Float64}())
    @test PolyJuMP.Bridges.Constraint.ToPolynomialBridge{Float64} in list
    @test PolyJuMP.Bridges.Objective.ToPolynomialBridge{Float64} in list
end

function test_optimizer_unconstrained()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @objective(model, Min, a^4 - 2a^2 + 1)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ 0 atol = 1e-6
    # The two minimizers `±1` are the atoms of the moment matrix
    @test result_count(model) == 2
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test sort([value(a; result = i) for i in 1:2]) ≈ [-1, 1] rtol = 1e-3
    @test objective_value(model) ≈ 0 atol = 1e-6
end

function test_optimizer_multiplier_maxdegree()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @objective(model, Min, a)
    @constraint(model, con, 1 - a^2 >= 0)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ -1 atol = 1e-6
    # `t` and the constant multiplier of `con`
    @test num_variables(unsafe_backend(model).relaxation) == 2
    # The minimizer `-1` is the atom of the moment matrix
    @test result_count(model) == 1
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(a) ≈ -1 rtol = 1e-3
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), con, 2)
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), con) == 2
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ -1 atol = 1e-6
    # `t` and the quadratic multiplier of `con`
    @test num_variables(unsafe_backend(model).relaxation) == 4
    @test result_count(model) == 1
    @test value(a) ≈ -1 rtol = 1e-3
end

function test_optimizer_equality()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @objective(model, Min, a)
    @constraint(model, a^2 == 1)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ -1 atol = 1e-6
    # The minimizer `-1` is the atom of the moment matrix
    @test result_count(model) == 1
    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(a) ≈ -1 rtol = 1e-3
    @test objective_value(model) ≈ -1 rtol = 1e-3
end

function test_algebraic_set_without_groebner_basis()
    @polyvar x y z
    # A single equation is a Gröbner basis so we add a redundant one
    K = @set x^2 + y^2 + z^2 == 1 && 2x^2 + 2y^2 + 2z^2 == 2
    model = SOSModel(SOLVER)
    @variable(model, t)
    @objective(model, Max, t)
    # `x + 1 = ((x + 1)^2 + y^2 + z^2) / 2` on the sphere
    @constraint(model, x - t in SOSCone(), domain = K, maxdegree = 2)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(t) ≈ -1 rtol = 1e-6
    # The equations are obtained from the Macaulay matrix instead
    @test !SemialgebraicSets.is_gröbner_basis(ideal(K))
    return
end

function test_algebraic_set_sparse_equations()
    # The dense equations obtained without the reduced row echelon form
    # gave `NUMERICAL_ERROR` with Clarabel
    @polyvar x y z
    # A single equation is a Gröbner basis so we add a redundant one
    K = @set x^2 + y^2 + z^2 == 1 && 2x^2 + 2y^2 + 2z^2 == 2
    p = x^4 * y^2 + x^2 * y^4 + z^6 - 3x^2 * y^2 * z^2 + 1
    model = SOSModel(SOLVER)
    @constraint(model, p in SOSCone(), domain = K)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test !SemialgebraicSets.is_gröbner_basis(ideal(K))
    return
end

# Maximize `t` such that `p - t` is SOS on `domain`
function _max_lower_bound(p, domain, maxdegree...)
    model = SOSModel(SOLVER)
    @variable(model, t)
    @objective(model, Max, t)
    if isempty(maxdegree)
        @constraint(model, p - t in SOSCone(), domain = domain)
    else
        @constraint(
            model,
            p - t in SOSCone(),
            domain = domain,
            maxdegree = only(maxdegree),
        )
    end
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    return value(t)
end

# With a smaller degree, the problems are weakly infeasible so we do not test
# them as the solver may not detect it.
function test_truncated_algebraic_set()
    @polyvar x y
    # `x - y = y * (x * y - 1) - x * (y^2 - 1)` vanishes on `V` but it is only
    # in the truncation of the ideal at degree 3
    V = @set x * y == 1 && y^2 == 1
    # The ideal is truncated at `maxdegree` and not at the degree of the
    # zero basis which is 2
    @test _max_lower_bound(x - y, V, 3) ≈ 0 atol = 1e-6
    # The degree of the truncation is independent of `maxdegree` which is 2
    # by default
    T = SemialgebraicSets.truncated(V, 3)
    @test _max_lower_bound(x - y, T) ≈ 0 atol = 1e-6
    @test _max_lower_bound(x - y, T, 2) ≈ 0 atol = 1e-6
    model = SOSModel(SOLVER)
    @variable(model, t)
    for domain in [V, T]
        @test_throws ArgumentError @constraint(
            model,
            x - y - t in SOSCone(),
            domain = domain,
            newton_of_remainder = true,
        )
    end
    @test !SemialgebraicSets.is_gröbner_basis(ideal(V))
    return
end

function test_truncated_algebraic_set_with_inequalities()
    @polyvar x y
    # On `K`, `x - 1 = (x - y) + (y - 1)` where `x - y` is only in the
    # truncation of the ideal at degree 3
    V = @set x * y == 1 && y^2 == 1
    K = V ∩ @set(y - 1 >= 0)
    @test _max_lower_bound(x, K, 3) ≈ 1 atol = 1e-6
    # The multipliers of the inequalities have degree at most 2 while the
    # multipliers of the equalities have degree at most 3
    K = SemialgebraicSets.truncated(V, 3) ∩ @set(y - 1 >= 0)
    @test _max_lower_bound(x, K) ≈ 1 atol = 1e-6
    @test !SemialgebraicSets.is_gröbner_basis(ideal(V))
    return
end

function test_truncated_newton()
    @polyvar x y z
    # The equation is redundant so that it is not a Gröbner basis
    V = @set x^2 + y^2 + z^2 == 1 && 2x^2 + 2y^2 + 2z^2 == 2
    for (K, bound) in [(V, -1.1720537), (V ∩ @set(x >= 0), -1)]
        sizes = map([true, false]) do newton
            model = SOSModel(SOLVER)
            @variable(model, t)
            @objective(model, Max, t)
            if newton
                c = @constraint(
                    model,
                    x^3 + y - t in SOSCone(),
                    domain = K,
                    maxdegree = 4,
                )
            else
                c = @constraint(
                    model,
                    x^3 + y - t in SOSCone(),
                    domain = K,
                    maxdegree = 4,
                    newton_polytope = nothing,
                )
            end
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(t) ≈ bound rtol = 1e-5
            return length(gram_matrix(c).basis)
        end
        # `x^2` is the leading monomial of `x^2 + y^2 + z^2 - 1`
        @test sizes == [9, 10]
    end
    @test !SemialgebraicSets.is_gröbner_basis(ideal(V))
    return
end

function _sparse_lower_bound(p, domain, sparsity, maxdegree = 2)
    model = SOSModel(SOLVER)
    @variable(model, t)
    @objective(model, Max, t)
    c = @constraint(
        model,
        p - t in SOSCone(),
        domain = domain,
        sparsity = sparsity,
        maxdegree = maxdegree,
    )
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    return value(t), c
end

function test_sparsity_with_equalities()
    @polyvar x y
    # On the domain, `y + 1/4 = (x + 1/2)^2` so the bound is `-1/4`.
    # `y` is invariant under `x ↦ -x` but `y - x^2 - x` is not so the Gram
    # matrix should not be block-diagonalized with this sign symmetry.
    V = @set y == x^2 + x
    # The equation is redundant so that it is not a Gröbner basis
    W = @set y == x^2 + x && 2y == 2x^2 + 2x
    for domain in [V, W, SemialgebraicSets.truncated(W, 2)]
        for sparsity in
            [Sparsity.SignSymmetry(), Sparsity.Monomial(), Sparsity.Variable()]
            bound, _ = _sparse_lower_bound(y, domain, sparsity)
            @test bound ≈ -1 / 4 rtol = 1e-5
        end
    end
    @test !SemialgebraicSets.is_gröbner_basis(ideal(W))
    return
end

function test_equality_multiplier_term_sparsity()
    @polyvar x y z
    # Multiplying the equality by z² gives the certificate (xz + z/2)².
    # Its cross term x*z² is absent from both the objective and the equality.
    V = SemialgebraicSets.truncated(@set(y == x^2 + x), 4)
    for (K, pattern) in [
        (V, Sparsity.Monomial()),
        (
            V ∩ @set(1 - z^2 >= 0),
            Sparsity.Monomial(ChordalCompletion(), 1, true),
        ),
    ]
        bound, c = _sparse_lower_bound(z^2 * (y + 0.25), K, pattern, 4)
        @test bound ≈ 0 atol = 1e-6
        @test any(gram_matrix(c).blocks) do block
            monos = MB.keys_as_monomials(block.basis)
            return z in monos && x * z in monos
        end
    end
    return
end

function test_scalar_localizing_constraints()
    @polyvar x y
    circle = SemialgebraicSets.truncated(@set(x^2 + y^2 == 1), 2)
    exterior = @set x^2 + y^2 >= 1
    for domain in (circle, exterior)
        bound, c = _sparse_lower_bound(x^2 + 2y^2, domain, Sparsity.Variable())
        @test bound ≈ 1 atol = 1e-6
        blocks = gram_matrix(c).blocks
        @test Set(MB.keys_as_monomials(b.basis) for b in blocks) ==
              Set(monomial_vector.([[1, x], [1, y]]))
    end
    return
end

function test_groebner_standard_monomials()
    @polyvar x
    # `x^100 = 1` on the domain so the lower bound of `x^100 + x` is `0` and
    # the default `maxdegree` is 100 but only the standard monomials `1` and
    # `x` are needed in the Gram basis
    model = SOSModel(SOLVER)
    @variable(model, t)
    @objective(model, Max, t)
    c = @constraint(model, x^100 + x - t in SOSCone(), domain = @set(x^2 == 1))
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(t) ≈ 0 atol = 1e-6
    @test MB.keys_as_monomials(gram_matrix(c).basis) == [1, x]
    return
end

# `min ∫_{-1}^{1} p` over the polynomials `p = a + b * t + c * t^2` that are
# nonnegative with `∫_{-1}^{1} p^2 = 1`. Substituting `p = λ(t - r)^2` and
# setting the derivative in `u = r^2` to zero gives `(5u - 1)^2 = 0` hence the
# optimal value `2√2/3` attained at `r = 1/√5`.
const _OPT = 2sqrt(2) / 3

function _bound_model(matrix_degree, equality_degree)
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    # the Gram matrix of `p` in the basis `(1, t)`
    matrix = @constraint(model, [a b/2; b/2 c] in PSDCone())
    equality = @constraint(
        model,
        2a^2 + (2 // 3) * (b^2 + 2a * c) + (2 // 5) * c^2 == 1
    )
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), matrix, matrix_degree)
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), equality, equality_degree)
    @objective(model, Min, 2a + (2 // 3) * c)
    optimize!(model)
    return model, a, b, c, matrix
end

function test_matrix_inequality_degree_0()
    model, _, _, _, matrix = _bound_model(0, 0)
    @test termination_status(model) == MOI.OPTIMAL
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), matrix) == 0
    # A constant multiplier only gives the trivial bound implied by `p ≥ 0`
    @test objective_bound(model) ≈ 0 atol = 1e-6
    @test result_count(model) == 0
end

function test_matrix_inequality_degree_2()
    model, a, b, c, _ = _bound_model(2, 2)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ _OPT rtol = 1e-5
    # `b` only appears squared in the constraints and not in the objective so
    # both `±b` are optimal and are the two atoms of the moment matrix
    @test result_count(model) == 2
    λ = 1 / sqrt(1.28)
    for i in 1:2
        @test primal_status(model; result = i) == MOI.FEASIBLE_POINT
        @test objective_value(model; result = i) ≈ _OPT rtol = 1e-4
        @test value(a; result = i) ≈ λ / 5 rtol = 1e-3
        @test abs(value(b; result = i)) ≈ 2λ / sqrt(5) rtol = 1e-3
        @test value(c; result = i) ≈ λ rtol = 1e-3
    end
    @test value(b; result = 1) ≈ -value(b; result = 2) rtol = 1e-3
end

# `Symmetric` gives `MOI.PositiveSemidefiniteConeTriangle` instead of
# `MOI.PositiveSemidefiniteConeSquare`
function test_matrix_inequality_triangle()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    matrix = @constraint(model, Symmetric([a b/2; b/2 c]) in PSDCone())
    equality = @constraint(
        model,
        2a^2 + (2 // 3) * (b^2 + 2a * c) + (2 // 5) * c^2 == 1
    )
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), matrix, 2)
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), equality, 2)
    @objective(model, Min, 2a + (2 // 3) * c)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ _OPT rtol = 1e-5
end

# The square cone also constrains the matrix to be symmetric so the entries
# that are not syntactically equal give an equality constraint
function test_matrix_inequality_square_symmetry()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    @variable(model, d)
    matrix = @constraint(model, [a b/2; d c] in PSDCone())
    equality = @constraint(
        model,
        2a^2 + (2 // 3) * (b^2 + 2a * c) + (2 // 5) * c^2 == 1
    )
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), matrix, 2)
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), equality, 2)
    @objective(model, Min, 2a + (2 // 3) * c)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ _OPT rtol = 1e-5
end

# The Sum-of-Squares bridges reformulate the `SOSCone` constraint into a
# matrix inequality on its Gram matrix
function test_matrix_inequality_from_sos_constraint()
    DynamicPolynomials.@polyvar t
    function integral(p)
        q = antidifferentiate(p, t)
        return coefficients(subs(q, t => 1))[1] -
               coefficients(subs(q, t => -1))[1]
    end
    model = SOSModel(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    p = a + b * t + c * t^2
    @constraint(model, p >= 0)
    @constraint(model, integral(p^2) == 1)
    @objective(model, Min, integral(p))
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    # a valid lower bound, the default degrees are too low for it to be tight
    @test objective_bound(model) <= _OPT + 1e-6
end

# A bilinear matrix inequality: `min x + y` such that `x^2 + y^2 == 4` and
# `|x * y| ≤ 1`. With `s = x + y` and `p = x * y`, `s^2 = 4 + 2p ≤ 6` so the
# optimal value is `-√6`, attained with `p = 1`.
function test_bilinear_matrix_inequality()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, x)
    @variable(model, y)
    matrix = @constraint(model, [1 x*y; x*y 1] in PSDCone())
    equality = @constraint(model, x^2 + y^2 == 4)
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), matrix, 2)
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), equality, 2)
    @objective(model, Min, x + y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ -sqrt(6) rtol = 1e-5
    # `x` and `y` can be swapped so there are two minimizers
    @test result_count(model) == 2
    for i in 1:2
        @test primal_status(model; result = i) == MOI.FEASIBLE_POINT
        @test value(x; result = i) + value(y; result = i) ≈ -sqrt(6) rtol = 1e-4
        @test value(x; result = i) * value(y; result = i) ≈ 1 rtol = 1e-4
    end
end

# Entries of degree at least three give a `MOI.VectorNonlinearFunction`:
# `|x^3| ≤ 1` so the optimal value is `-1` for `min x` and `1` for `max x`
function _cubic_model(sense)
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, x)
    matrix = @constraint(model, [1 x^3; x^3 1] in PSDCone())
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), matrix, 2)
    if sense != MOI.FEASIBILITY_SENSE
        @objective(model, sense, x)
    end
    optimize!(model)
    return model, x
end

function test_cubic_matrix_inequality()
    for (sense, opt) in [(MIN_SENSE, -1), (MAX_SENSE, 1)]
        model, x = _cubic_model(sense)
        @test termination_status(model) == MOI.OPTIMAL
        @test objective_bound(model) ≈ opt rtol = 1e-5
        @test result_count(model) == 1
        @test primal_status(model) == MOI.FEASIBLE_POINT
        @test value(x) ≈ opt rtol = 1e-4
    end
end

function test_cubic_matrix_inequality_feasibility()
    model, _ = _cubic_model(MOI.FEASIBILITY_SENSE)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ 0 atol = 1e-6
end

function test_matrix_inequality_empty()
    optimizer = SumOfSquares.Optimizer(SOLVER)
    @test MOI.is_empty(optimizer)
    x = MOI.add_variables(optimizer, 3)
    ci = MOI.add_constraint(
        optimizer,
        MOI.VectorOfVariables(x),
        MOI.PositiveSemidefiniteConeTriangle(2),
    )
    @test MOI.is_valid(optimizer, ci)
    @test !MOI.is_empty(optimizer)
    @test MOI.supports(optimizer, PolyJuMP.MultiplierMaxdegree(), typeof(ci))
    @test isnothing(MOI.get(optimizer, PolyJuMP.MultiplierMaxdegree(), ci))
    MOI.set(optimizer, PolyJuMP.MultiplierMaxdegree(), ci, 4)
    @test MOI.get(optimizer, PolyJuMP.MultiplierMaxdegree(), ci) == 4
    MOI.empty!(optimizer)
    @test MOI.is_empty(optimizer)
    @test !MOI.is_valid(optimizer, ci)
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end # module

TestOptimizer.runtests()
