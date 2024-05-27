using Test
using SumOfSquares
using DynamicPolynomials

# Like `term_test` but with a variable `y` fixed to 1.
function term_fixed_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
)
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
    @test objective_value(model) ≈ 1.0 atol = atol rtol = rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 1.0 atol = atol rtol = rtol

    test_constraint_primal(cref, 0.0)

    p = gram_matrix(cref)
    @test value_matrix(p) ≈ zeros(1, 1) atol = atol rtol = rtol
    @test p.basis.monomials == [x]

    @test dual_status(model) == MOI.FEASIBLE_POINT
    for (m, μ) in [(x^2, moments(cref))]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 1
        @test moment_value(moments(μ)[1]) ≈ 1.0 atol = atol rtol = rtol
        @test moments(μ)[1].polynomial.monomial == m
    end

    ν = moment_matrix(cref)
    @test value_matrix(ν) ≈ ones(1, 1) atol = atol rtol = rtol
    @test ν.basis.monomials == [x]

    N = SumOfSquares.Certificate.NewtonFilter{
        SumOfSquares.Certificate.NewtonDegreeBounds{Tuple{}},
    }
    S = SumOfSquares.SOSPolynomialSet{
        typeof(set),
        SubBasis{Monomial,monomial_type(x),monomial_vector_type(x)},
        SumOfSquares.Certificate.Remainder{
            SumOfSquares.Certificate.Newton{
                typeof(cone),
                FullBasis{Monomial,monomial_type(x)},
                N,
            },
        },
    }
    @test list_of_constraint_types(model) == [(Vector{JuMP.AffExpr}, S)]
    return test_delete_bridge(
        model,
        cref,
        1,
        (
            (MOI.VectorOfVariables, MOI.Nonnegatives, 0),
            (
                MOI.VectorAffineFunction{Float64},
                SumOfSquares.PolyJuMP.ZeroPolynomialSet{
                    typeof(set),
                    SubBasis{Monomial,monomial_type(x),monomial_vector_type(x)},
                },
                0,
            ),
        ),
    )
end
function sos_term_fixed_test(optimizer, config)
    return term_fixed_test(optimizer, config, SOSCone())
end
sd_tests["sos_term_fixed"] = sos_term_fixed_test
function sdsos_term_fixed_test(optimizer, config)
    return term_fixed_test(optimizer, config, SDSOSCone())
end
soc_tests["sdsos_term_fixed"] = sdsos_term_fixed_test
function dsos_term_fixed_test(optimizer, config)
    return term_fixed_test(optimizer, config, DSOSCone())
end
linear_tests["dsos_term_fixed"] = dsos_term_fixed_test
