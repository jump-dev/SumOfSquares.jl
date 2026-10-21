module TestOptimizer

using Test

using LinearAlgebra
import MathOptInterface as MOI
using JuMP
using DynamicPolynomials
using SumOfSquares

import Clarabel
const SOLVER =
    optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

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
