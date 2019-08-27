using Test
using SumOfSquares
using DynamicPolynomials

function quadratic_test(
    optimizer, config::MOIT.TestConfig,
    cone::SumOfSquares.PolyJuMP.PolynomialSet, bivariate::Bool)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x
    if bivariate
        @polyvar y
        poly = x^2 + α*x*y + y^2
        cert_monos = [x, y]
        monos = [x^2, x*y, y^2]
    else
        poly = x^2 + α*x + 1
        cert_monos = [x, 1]
        monos = [x^2, x, 1]
    end
    cref = @constraint(model, poly in cone)

    # See https://github.com/JuliaOpt/MathOptInterface.jl/issues/676
    @objective(model, Max, α + 1)
    optimize!(model)

    @test certificate_monomials(cref) == cert_monos

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 3.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 2.0 atol=atol rtol=rtol

    @test_throws SumOfSquares.ValueNotSupported value(cref)
    p = gram_matrix(cref)
    @test getmat(p) ≈ ones(2, 2) atol=atol rtol=rtol
    @test p.x == cert_monos

    a = moment_value.(moments(dual(cref)))
    @test a[2] ≈ -1.0 atol=atol rtol=rtol
    @test a[1] + a[3] ≈ 2.0 atol=atol rtol=rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    for μ in [dual(cref), moments(cref)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 3
        @test a ≈ moment_value.(moments(μ)) atol=atol rtol=rtol
        @test monomial.(moments(μ)) == monos
    end

    ν = moment_matrix(cref)
    @test getmat(ν) ≈ [a[1] a[2]
                       a[2] a[3]] atol=atol rtol=rtol
    @test ν.x == cert_monos

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
sos_univariate_quadratic_test(optimizer, config)   = quadratic_test(optimizer, config, SOSCone(), false)
sd_tests["sos_univariate_quadratic"] = sos_univariate_quadratic_test
sos_bivariate_quadratic_test(optimizer, config)    = quadratic_test(optimizer, config, SOSCone(), true)
sd_tests["sos_bivariate_quadratic"] = sos_bivariate_quadratic_test
sdsos_univariate_quadratic_test(optimizer, config) = quadratic_test(optimizer, config, SDSOSCone(), false)
soc_tests["sdsos_univariate_quadratic"] = sdsos_univariate_quadratic_test
sdsos_bivariate_quadratic_test(optimizer, config)  = quadratic_test(optimizer, config, SDSOSCone(), true)
soc_tests["sdsos_bivariate_quadratic"] = sdsos_bivariate_quadratic_test
dsos_univariate_quadratic_test(optimizer, config)  = quadratic_test(optimizer, config, DSOSCone(), false)
linear_tests["dsos_univariate_quadratic"] = dsos_univariate_quadratic_test
dsos_bivariate_quadratic_test(optimizer, config)   = quadratic_test(optimizer, config, DSOSCone(), true)
linear_tests["dsos_bivariate_quadratic"] = dsos_bivariate_quadratic_test
