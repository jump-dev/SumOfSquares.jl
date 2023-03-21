using Test
using SumOfSquares
using DynamicPolynomials

function term_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x
    cref = @constraint(model, α * x^2 in cone)

    # See https://github.com/jump-dev/MathOptInterface.jl/issues/676
    @objective(model, Min, α + 1)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 1.0 atol = atol rtol = rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 0.0 atol = atol rtol = rtol

    test_constraint_primal(cref, 0.0)

    p = gram_matrix(cref)
    @test getmat(p) ≈ zeros(1, 1) atol = atol rtol = rtol
    @test p.basis.monomials == [x]

    @test dual_status(model) == MOI.FEASIBLE_POINT
    for μ in [dual(cref), moments(cref)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 1
        @test moment_value(moments(μ)[1]) ≈ 1.0 atol = atol rtol = rtol
        @test monomial(moments(μ)[1]) == x^2
    end

    ν = moment_matrix(cref)
    @test getmat(ν) ≈ ones(1, 1) atol = atol rtol = rtol
    @test ν.basis.monomials == [x]

    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace,
        Monomial{true},
        MonomialVector{true},
        SumOfSquares.Certificate.Newton{typeof(cone),MonomialBasis,Tuple{}},
    }
    @test list_of_constraint_types(model) == [(Vector{VariableRef}, S)]
    return test_delete_bridge(
        model,
        cref,
        1,
        (
            (MOI.VectorOfVariables, MOI.Nonnegatives, 0),
            (
                MOI.VectorAffineFunction{Float64},
                SumOfSquares.PolyJuMP.ZeroPolynomialSet{
                    SumOfSquares.FullSpace,
                    MonomialBasis,
                    Monomial{true},
                    MonomialVector{true},
                },
                0,
            ),
        ),
    )
end
sos_term_test(optimizer, config) = term_test(optimizer, config, SOSCone())
sd_tests["sos_term"] = sos_term_test
sdsos_term_test(optimizer, config) = term_test(optimizer, config, SDSOSCone())
soc_tests["sdsos_term"] = sdsos_term_test
dsos_term_test(optimizer, config) = term_test(optimizer, config, DSOSCone())
linear_tests["dsos_term"] = dsos_term_test
