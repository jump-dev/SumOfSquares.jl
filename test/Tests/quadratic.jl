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

    @test certificate_monomials(cref) == [x, 1]

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 3.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 2.0 atol=atol rtol=rtol

    @test_throws SumOfSquares.ValueNotSupported value(cref)
    p = gram_matrix(cref)
    @test getmat(p) ≈ ones(2, 2) atol=atol rtol=rtol
    @test p.x == [x, 1]

    a = moment_value.(moments(dual(cref)))
    @test a[2] ≈ -1.0 atol=atol rtol=rtol
    @test a[1] + a[3] ≈ 2.0 atol=atol rtol=rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    for μ in [dual(cref), moments(cref)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 3
        @test a ≈ moment_value.(moments(μ)) atol=atol rtol=rtol
        @test monomial.(moments(μ)) == [x^2, x, 1]
    end

    ν = moment_matrix(cref)
    @test getmat(ν) ≈ [a[1] a[2]
                       a[2] a[3]] atol=atol rtol=rtol
    @test ν.x == [x, 1]

    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace, Monomial{true}, MonomialVector{true}, SumOfSquares.Certificate.Remainder{typeof(cone), SumOfSquares.MonomialBasis, Tuple{}}
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
sd_tests["sos_bivariate_quadratic"] = sos_bivariate_quadratic_test
sdsos_bivariate_quadratic_test(optimizer, config) = bivariate_quadratic_test(optimizer, config, SDSOSCone())
soc_tests["sdsos_bivariate_quadratic"] = sdsos_bivariate_quadratic_test
dsos_bivariate_quadratic_test(optimizer, config)  = bivariate_quadratic_test(optimizer, config, DSOSCone())
linear_tests["dsos_bivariate_quadratic"] = dsos_bivariate_quadratic_test
