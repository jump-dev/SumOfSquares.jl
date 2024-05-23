using Test
import MultivariateBases
using DynamicPolynomials

function quadratic_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
    basis,
    bivariate::Bool,
)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @variable(model, α)

    @polyvar x
    if bivariate
        @polyvar y
        poly = x^2 + α * x * y + y^2
        cert_monos = [y, x]
        monos = [y^2, x * y, x^2]
    else
        poly = x^2 + α * x + 1
        cert_monos = [1, x]
        monos = [1, x, x^2]
    end
    cref = @constraint(model, poly in cone, basis = basis)

    @objective(model, Max, α)
    optimize!(model)

    if basis == ChebyshevBasis
        err = ErrorException(
            "`certificate_monomials` is not supported with `$(basis{typeof(x + 1.0)})`, use `certificate_basis` instead.",
        )
        @test_throws err certificate_monomials(cref)
    else
        @test certificate_monomials(cref) == cert_monos
    end

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 2.0 atol = atol rtol = rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 2.0 atol = atol rtol = rtol

    test_constraint_primal(cref, value(poly))

    p = gram_matrix(cref)
    @test value_matrix(p) ≈ ones(2, 2) atol = atol rtol = rtol
    if basis == ChebyshevBasis
        @test p.basis.polynomials == cert_monos
    else
        @test p.basis.monomials == cert_monos
    end

    a = moment_value.(moments(dual(cref)))
    @test a[2] ≈ -1.0 atol = atol rtol = rtol
    @test a[1] + a[3] ≈ 2.0 atol = atol rtol = rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    for μ in [dual(cref), moments(cref)]
        @test μ isa AbstractMeasure{Float64}
        @test length(moments(μ)) == 3
        @test a ≈ moment_value.(moments(μ)) atol = atol rtol = rtol
        @test monomial.(moments(μ)) == monos
    end

    ν = moment_matrix(cref)
    @test value_matrix(ν) ≈ [
        a[1] a[2]
        a[2] a[3]
    ] atol = atol rtol = rtol
    if basis == ChebyshevBasis
        @test p.basis.polynomials == cert_monos
    else
        @test ν.basis.monomials == cert_monos
    end

    N = SumOfSquares.Certificate.NewtonFilter{
        SumOfSquares.Certificate.NewtonDegreeBounds{Tuple{}},
    }
    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace,
        MB.SubBasis{MB.Monomial,monomial_type(x),monomial_vector_type(x)},
        SumOfSquares.Certificate.Newton{typeof(cone),basis,N},
    }
    @test list_of_constraint_types(model) == [(Vector{AffExpr}, S)]
    return test_delete_bridge(
        model,
        cref,
        1,
        (
            (MOI.VectorOfVariables, MOI.Nonnegatives, 0),
            (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}, 0),
            (
                MOI.VectorAffineFunction{Float64},
                SumOfSquares.PolyJuMP.ZeroPolynomialSet{
                    SumOfSquares.FullSpace,
                    basis,
                    monomial_type(x),
                    monomial_vector_type(x),
                },
                0,
            ),
        ),
    )
end
function sos_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), MonomialBasis, false)
end
sd_tests["sos_univariate_quadratic"] = sos_univariate_quadratic_test
function sos_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), MonomialBasis, true)
end
sd_tests["sos_bivariate_quadratic"] = sos_bivariate_quadratic_test
function sos_scaled_univariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        SOSCone(),
        ScaledMonomialBasis,
        false,
    )
end
sd_tests["sos_scaled_univariate_quadratic"] =
    sos_scaled_univariate_quadratic_test
function sos_scaled_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        SOSCone(),
        ScaledMonomialBasis,
        true,
    )
end
sd_tests["sos_scaled_bivariate_quadratic"] = sos_scaled_bivariate_quadratic_test
function sos_cheby_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), ChebyshevBasis, false)
end
sd_tests["sos_cheby_univariate_quadratic"] = sos_cheby_univariate_quadratic_test
function sos_cheby_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), ChebyshevBasis, true)
end
sd_tests["sos_scaled_bivariate_quadratic"] = sos_scaled_bivariate_quadratic_test

function sdsos_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), MonomialBasis, false)
end
soc_tests["sdsos_univariate_quadratic"] = sdsos_univariate_quadratic_test
function sdsos_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), MonomialBasis, true)
end
soc_tests["sdsos_bivariate_quadratic"] = sdsos_bivariate_quadratic_test
function sdsos_scaled_univariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        SDSOSCone(),
        ScaledMonomialBasis,
        false,
    )
end
soc_tests["sdsos_scaled_univariate_quadratic"] =
    sdsos_scaled_univariate_quadratic_test
function sdsos_scaled_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        SDSOSCone(),
        ScaledMonomialBasis,
        true,
    )
end
soc_tests["sdsos_scaled_bivariate_quadratic"] =
    sdsos_scaled_bivariate_quadratic_test
function sdsos_cheby_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), ChebyshevBasis, false)
end
soc_tests["sdsos_cheby_univariate_quadratic"] =
    sdsos_cheby_univariate_quadratic_test
function sdsos_cheby_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), ChebyshevBasis, true)
end
soc_tests["sdsos_scaled_bivariate_quadratic"] =
    sdsos_scaled_bivariate_quadratic_test

function dsos_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), MonomialBasis, false)
end
linear_tests["dsos_univariate_quadratic"] = dsos_univariate_quadratic_test
function dsos_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), MonomialBasis, true)
end
linear_tests["dsos_bivariate_quadratic"] = dsos_bivariate_quadratic_test
function dsos_scaled_univariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        DSOSCone(),
        ScaledMonomialBasis,
        false,
    )
end
linear_tests["dsos_scaled_univariate_quadratic"] =
    dsos_scaled_univariate_quadratic_test
function dsos_scaled_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        DSOSCone(),
        ScaledMonomialBasis,
        true,
    )
end
linear_tests["dsos_scaled_bivariate_quadratic"] =
    dsos_scaled_bivariate_quadratic_test
function dsos_cheby_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), ChebyshevBasis, false)
end
linear_tests["dsos_cheby_univariate_quadratic"] =
    dsos_cheby_univariate_quadratic_test
function dsos_cheby_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), ChebyshevBasis, true)
end
linear_tests["dsos_cheby_bivariate_quadratic"] =
    dsos_cheby_bivariate_quadratic_test
