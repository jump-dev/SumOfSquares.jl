using Test
using SumOfSquares
using DynamicPolynomials

# Like `term_test` but with a variable `y` fixed to 1.
function term_fixed_test(
    optimizer, config::MOIT.TestConfig, cone::SumOfSquares.PolyJuMP.PolynomialSet)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x y
    set = @set y == 1
    cref = @constraint(model, (α - 1) * x^2 * y in cone, domain = set)

    @objective(model, Min, α)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 1.0 atol=atol rtol=rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 1.0 atol=atol rtol=rtol

    @test_throws SumOfSquares.ValueNotSupported value(cref)
    p = gram_matrix(cref)
    @test getmat(p) ≈ zeros(1, 1) atol=atol rtol=rtol
    @test p.x == [x]

    @test dual_status(model) == MOI.FEASIBLE_POINT
    for (m, μ) in [(x^2 * y, dual(cref)), (x^2, moments(cref))]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 1
        @test moment_value(moments(μ)[1]) ≈ 1.0 atol=atol rtol=rtol
        @test monomial(moments(μ)[1]) == m
    end

    ν = moment_matrix(cref)
    @test getmat(ν) ≈ ones(1, 1) atol=atol rtol=rtol
    @test ν.x == [x]

    S = SumOfSquares.SOSPolynomialSet{
        typeof(set), Monomial{true}, MonomialVector{true}, SumOfSquares.Certificate.Remainder{typeof(cone), SumOfSquares.MonomialBasis, Tuple{}}
    }
    @test list_of_constraint_types(model) == [(Vector{JuMP.AffExpr}, S)]
    test_delete_bridge(
        model, cref, 1,
        ((MOI.VectorOfVariables, MOI.Nonnegatives, 0),
         (MOI.VectorAffineFunction{Float64},
          SumOfSquares.PolyJuMP.ZeroPolynomialSet{
              typeof(set), SumOfSquares.MonomialBasis,
              Monomial{true}, MonomialVector{true}},
          0)))
end
sos_term_fixed_test(optimizer, config)   = term_fixed_test(optimizer, config, SOSCone())
sd_tests["sos_term_fixed"] = sos_term_fixed_test
sdsos_term_fixed_test(optimizer, config) = term_fixed_test(optimizer, config, SDSOSCone())
soc_tests["sdsos_term_fixed"] = sdsos_term_fixed_test
dsos_term_fixed_test(optimizer, config)  = term_fixed_test(optimizer, config, DSOSCone())
linear_tests["dsos_term_fixed"] = dsos_term_fixed_test
