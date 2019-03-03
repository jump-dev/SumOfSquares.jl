using Test
using SumOfSquares
using DynamicPolynomials

function bivariate_quadratic_test(optimizer,
                                  config::MOIT.TestConfig,
                                  cone::SumOfSquares.PolyJuMP.PolynomialSet)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x
    cref = @constraint(model, x^2 + α*x + 1 in cone)

    # See https://github.com/JuliaOpt/MathOptInterface.jl/issues/676
    @objective(model, Max, α + 1)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 3.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 2.0 atol=atol rtol=rtol

    @test_throws SumOfSquares.ValueNotSupported value(cref)
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

    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace, typeof(cone), SumOfSquares.MonomialBasis,
        Monomial{true},MonomialVector{true},Tuple{}
    }
    @test list_of_constraint_types(model) == [(Vector{AffExpr}, S)]
    test_delete_bridge(
        model, cref, 1,
        ((MOI.VectorOfVariables, MOI.Nonnegatives, 0),
         (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}, 0),
         (MOI.VectorAffineFunction{Float64},
          SumOfSquares.PolyJuMP.ZeroPolynomialSet{
              SumOfSquares.FullSpace, SumOfSquares.MonomialBasis,
              Monomial{true}, MonomialVector{true}},
          0)))
end
sos_bivariate_quadratic_test(optimizer, config)   = bivariate_quadratic_test(optimizer, config, SOSCone())
sdsos_bivariate_quadratic_test(optimizer, config) = bivariate_quadratic_test(optimizer, config, SDSOSCone())
dsos_bivariate_quadratic_test(optimizer, config)  = bivariate_quadratic_test(optimizer, config, DSOSCone())
