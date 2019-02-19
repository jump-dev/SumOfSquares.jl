using Test
using SumOfSquares
using DynamicPolynomials

include("utilities.jl")

function bivariate_quadratic_test(optimizer,
                                  config::MOI.Test.TestConfig,
                                  cone::SumOfSquares.PolyJuMP.PolynomialSet)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x
    cref = @constraint(model, x^2 + α*x + 1 in cone)

    @objective(model, Max, α)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 2.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 2.0 atol=atol rtol=rtol

    p = gram_matrix(cref)
    @test getmat(p) ≈ ones(2, 2) atol=atol rtol=rtol
    @test p.x == [x, 1]

    @test dual_status(model) == MOI.FEASIBLE_POINT
    μ = dual(cref)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == 3
    a = moment_value.(moments(μ))
    @test a[2] ≈ -1.0  atol=atol rtol=rtol
    @test a[1] + a[3] ≈ 2.0  atol=atol rtol=rtol
    @test monomial.(moments(μ)) == [x^2, x, 1]

    ν = moment_matrix(cref)
    @test getmat(ν) ≈ [a[1] a[2]
                       a[2] a[3]] atol=atol rtol=rtol
    @test ν.x == [x, 1]
end
sos_bivariate_quadratic_test(optimizer, config)   = bivariate_quadratic_test(optimizer, config, SOSCone())
sdsos_bivariate_quadratic_test(optimizer, config) = bivariate_quadratic_test(optimizer, config, SDSOSCone())
dsos_bivariate_quadratic_test(optimizer, config)  = bivariate_quadratic_test(optimizer, config, DSOSCone())
