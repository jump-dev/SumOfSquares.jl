module TestOptimizer

using Test

import MathOptInterface as MOI
using JuMP
using SumOfSquares

import Clarabel
const SOLVER =
    optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)

function test_optimizer_attributes()
    optimizer = SumOfSquares.Optimizer(SOLVER)
    @test optimizer isa SumOfSquares.Optimizer{Float64}
    @test MOI.get(optimizer, MOI.SolverName()) == "SumOfSquares"
    @test MOI.get(optimizer, MOI.TerminationStatus()) ==
          MOI.OPTIMIZE_NOT_CALLED
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
    @test primal_status(model) == MOI.NO_SOLUTION
    @test result_count(model) == 0
    @test objective_bound(model) ≈ 0 atol = 1e-6
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
    MOI.set(model, PolyJuMP.MultiplierMaxdegree(), con, 2)
    @test MOI.get(model, PolyJuMP.MultiplierMaxdegree(), con) == 2
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ -1 atol = 1e-6
    # `t` and the quadratic multiplier of `con`
    @test num_variables(unsafe_backend(model).relaxation) == 4
end

function test_optimizer_equality()
    model = Model(() -> SumOfSquares.Optimizer(SOLVER))
    @variable(model, a)
    @objective(model, Min, a)
    @constraint(model, a^2 == 1)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_bound(model) ≈ -1 atol = 1e-6
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
