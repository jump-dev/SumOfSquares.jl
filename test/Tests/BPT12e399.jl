# Adapted from Example 3.99 of [BPT12].
#
# [BPT12] Blekherman, G.; Parrilo, P. & Thomas, R.
# Semidefinite Optimization and Convex Algebraic Geometry
# Society for Industrial and Applied Mathematics, 2012

using Test
using SumOfSquares
using DynamicPolynomials

function BPT12e399_test(optimizer, config::MOIT.TestConfig)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x y
    cref = @constraint(model, 10 - (x^2 + α*y) in SOSCone(),
                       domain = @set x^2 + y^2 == 1)

    @objective(model, Max, α)
    JuMP.optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 6.0 atol=atol rtol=rtol

    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.value(α) ≈ 6 atol=1e-6

    @test_throws SumOfSquares.ValueNotSupported value(cref)
    p = gram_matrix(cref)
    @test getmat(p) ≈ [1 -3; -3 9] atol=atol rtol=rtol
    @test p.x == [y, 1]

    @test dual_status(model) == MOI.FEASIBLE_POINT
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    @test moment_value(moments(μ)[1]) ≈ -8/3 atol=atol rtol=rtol
    @test monomial(moments(μ)[1]) == x^2
    @test moment_value(moments(μ)[2]) ≈ 1 atol=atol rtol=rtol
    @test monomial(moments(μ)[2]) == y
    @test moment_value(moments(μ)[3]) ≈ 1/3 atol=atol rtol=rtol
    @test monomial(moments(μ)[3]) == 1

    μ = moments(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    @test moment_value(moments(μ)[1]) ≈ 3 atol=atol rtol=rtol
    @test monomial(moments(μ)[1]) == y^2
    @test moment_value(moments(μ)[2]) ≈ 1 atol=atol rtol=rtol
    @test monomial(moments(μ)[2]) == y
    @test moment_value(moments(μ)[3]) ≈ 1/3 atol=atol rtol=rtol
    @test monomial(moments(μ)[3]) == 1

    @objective(model, Min, α)
    JuMP.optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ -6.0 atol=atol rtol=rtol

    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.value(α) ≈ -6 atol=1e-6

    @test_throws SumOfSquares.ValueNotSupported value(cref)
    p = gram_matrix(cref)
    @test getmat(p) ≈ [1 3; 3 9] atol=atol rtol=rtol
    @test p.x == [y, 1]

    @test dual_status(model) == MOI.FEASIBLE_POINT
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    @test moment_value(moments(μ)[1]) ≈ -8/3 atol=atol rtol=rtol
    @test monomial(moments(μ)[1]) == x^2
    @test moment_value(moments(μ)[2]) ≈ -1 atol=atol rtol=rtol
    @test monomial(moments(μ)[2]) == y
    @test moment_value(moments(μ)[3]) ≈ 1/3 atol=atol rtol=rtol
    @test monomial(moments(μ)[3]) == 1

    μ = moments(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    @test moment_value(moments(μ)[1]) ≈ 3 atol=atol rtol=rtol
    @test monomial(moments(μ)[1]) == y^2
    @test moment_value(moments(μ)[2]) ≈ -1 atol=atol rtol=rtol
    @test monomial(moments(μ)[2]) == y
    @test moment_value(moments(μ)[3]) ≈ 1/3 atol=atol rtol=rtol
    @test monomial(moments(μ)[3]) == 1
end

sd_tests["BPT12e399"] = BPT12e399_test
