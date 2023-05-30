# Adapted from Example 3.99 of [BPT12].
#
# [BPT12] Blekherman, G.; Parrilo, P. & Thomas, R.
# Semidefinite Optimization and Convex Algebraic Geometry
# Society for Industrial and Applied Mathematics, 2012

using Test
using SumOfSquares
using DynamicPolynomials

function BPT12e399_test(optimizer, config::MOI.Test.Config, remainder::Bool)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x y
    if remainder
        cref = @constraint(
            model,
            10 - (x^2 + α * y) in SOSCone(),
            newton_of_remainder = true,
            maxdegree = nothing,
            domain = @set x^2 + y^2 == 1
        )
    else
        cref = @constraint(
            model,
            10 - (x^2 + α * y) in SOSCone(),
            maxdegree = 2,
            domain = @set x^2 + y^2 == 1
        )
    end

    @objective(model, Max, α)
    JuMP.optimize!(model)

    α_value = remainder ? 6.0 : 10.0

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ α_value atol = atol rtol = rtol

    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.value(α) ≈ α_value atol = atol rtol = rtol

    test_constraint_primal(cref, 10 - (x^2 + α_value * y))

    p = gram_matrix(cref)
    if remainder
        @test value_matrix(p) ≈ [9 -3; -3 1] atol = atol rtol = rtol
        @test p.basis.monomials == [1, y]
    else
        @test value_matrix(p) ≈ [
             5 -5 0
            -5  5 0
             0  0 4
        ] atol = atol rtol = rtol
        @test p.basis.monomials == [1, y, x]
    end

    @test dual_status(model) == MOI.FEASIBLE_POINT
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    @test moment_value(moments(μ)[1]) ≈ (remainder ? 1 / 3 : 1.0) atol = atol rtol =
        rtol
    @test monomial(moments(μ)[1]) == 1
    @test moment_value(moments(μ)[2]) ≈ 1 atol = atol rtol = rtol
    @test monomial(moments(μ)[2]) == y
    @test moment_value(moments(μ)[3]) ≈ (remainder ? -8 / 3 : 0.0) atol = atol rtol =
        rtol
    @test monomial(moments(μ)[3]) == x^2

    μ = moments(cref)
    @test μ isa AbstractMeasure{Float64}
    if remainder
        @test length(moments(μ)) == 3
        @test moment_value(moments(μ)[1]) ≈ 1 / 3 atol = atol rtol = rtol
        @test monomial(moments(μ)[1]) == 1
        @test moment_value(moments(μ)[2]) ≈ 1 atol = atol rtol = rtol
        @test monomial(moments(μ)[2]) == y
        @test moment_value(moments(μ)[3]) ≈ 3 atol = atol rtol = rtol
        @test monomial(moments(μ)[3]) == y^2
    else
        @test length(moments(μ)) == 5
        @test moment_value(moments(μ)[1]) ≈ 1 atol = atol rtol = rtol
        @test monomial(moments(μ)[1]) == 1
        @test moment_value(moments(μ)[2]) ≈ 1 atol = atol rtol = rtol
        @test monomial(moments(μ)[2]) == y
        @test moment_value(moments(μ)[3]) ≈ 0 atol = atol rtol = rtol
        @test monomial(moments(μ)[3]) == x
        @test moment_value(moments(μ)[4]) ≈ 1 atol = atol rtol = rtol
        @test monomial(moments(μ)[4]) == y^2
        @test moment_value(moments(μ)[5]) ≈ 0 atol = atol rtol = rtol
        @test monomial(moments(μ)[5]) == x * y
    end

    @objective(model, Min, α)
    JuMP.optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ -α_value atol = atol rtol = rtol

    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.value(α) ≈ -α_value atol = atol rtol = rtol

    test_constraint_primal(cref, 10 - (x^2 - α_value * y))

    p = gram_matrix(cref)
    if remainder
        @test value_matrix(p) ≈ [9 3; 3 1] atol = atol rtol = rtol
        @test p.basis.monomials == [1, y]
    else
        @test value_matrix(p) ≈ [
            5 5 0
            5 5 0
            0 0 4
        ] atol = atol rtol = rtol
        @test p.basis.monomials == [1, y, x]
    end

    @test dual_status(model) == MOI.FEASIBLE_POINT
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    @test moment_value(moments(μ)[1]) ≈ (remainder ? 1 / 3 : 1.0) atol = atol rtol =
        rtol
    @test monomial(moments(μ)[1]) == 1
    @test moment_value(moments(μ)[2]) ≈ -1 atol = atol rtol = rtol
    @test monomial(moments(μ)[2]) == y
    @test moment_value(moments(μ)[3]) ≈ (remainder ? -8 / 3 : 0.0) atol = atol rtol =
        rtol
    @test monomial(moments(μ)[3]) == x^2

    μ = moments(cref)
    @test μ isa AbstractMeasure{Float64}
    if remainder
        @test length(moments(μ)) == 3
        @test moment_value(moments(μ)[1]) ≈ 1 / 3 atol = atol rtol = rtol
        @test monomial(moments(μ)[1]) == 1
        @test moment_value(moments(μ)[2]) ≈ -1 atol = atol rtol = rtol
        @test monomial(moments(μ)[2]) == y
        @test moment_value(moments(μ)[3]) ≈ 3 atol = atol rtol = rtol
        @test monomial(moments(μ)[3]) == y^2
    else
        @test length(moments(μ)) == 5
        @test moment_value(moments(μ)[1]) ≈ 1 atol = atol rtol = rtol
        @test monomial(moments(μ)[1]) == 1
        @test moment_value(moments(μ)[2]) ≈ -1 atol = atol rtol = rtol
        @test monomial(moments(μ)[2]) == y
        @test moment_value(moments(μ)[3]) ≈ 0 atol = atol rtol = rtol
        @test monomial(moments(μ)[3]) == x
        @test moment_value(moments(μ)[4]) ≈ 1 atol = atol rtol = rtol
        @test monomial(moments(μ)[4]) == y^2
        @test moment_value(moments(μ)[5]) ≈ 0 atol = atol rtol = rtol
        @test monomial(moments(μ)[5]) == x * y
    end

    return model
end

BPT12e399_rem_test(optimizer, config) = BPT12e399_test(optimizer, config, true)
sd_tests["BPT12e399_rem"] = BPT12e399_rem_test
function BPT12e399_maxdegree_test(optimizer, config)
    return BPT12e399_test(optimizer, config, false)
end
sd_tests["BPT12e399_maxdegree"] = BPT12e399_maxdegree_test
