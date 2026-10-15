module TestOptimizer

using Test

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
