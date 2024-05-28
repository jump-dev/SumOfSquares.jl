using Test
import MultivariateBases
using DynamicPolynomials

function _test_moments(μ, a, monos; atol, rtol)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == length(a)
    @test a ≈ moment_value.(moments(μ)) atol = atol rtol = rtol
    @test [m.polynomial.monomial for m in moments(μ)] == monos
end

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

    @test certificate_monomials(cref) == cert_monos

    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 2.0 atol = atol rtol = rtol

    @test primal_status(model) == MOI.FEASIBLE_POINT
    @test value(α) ≈ 2.0 atol = atol rtol = rtol

    test_constraint_primal(cref, value(poly))

    p = gram_matrix(cref)
    @test value_matrix(p) ≈ ones(2, 2) atol = atol rtol = rtol
    @test p.basis.monomials == cert_monos

    a = moment_value.(moments(dual(cref)))
    @test a[2] ≈ -1.0 atol = atol rtol = rtol
    @test a[1] + a[3] ≈ 2.0 atol = atol rtol = rtol

    @test dual_status(model) == MOI.FEASIBLE_POINT
    _test_moments(dual(cref), a, monos; atol, rtol)
    if basis === MB.Chebyshev && bivariate
        _test_moments(moments(cref), [0, 0, 2, -1, 0, 0], monomial_vector([1, x, y^2, x * y, x^2, x^3]); atol, rtol)
    else
        _test_moments(moments(cref), a, monos; atol, rtol)
    end

    ν = moment_matrix(cref)
    @test value_matrix(ν) ≈ [
        a[1] a[2]
        a[2] a[3]
    ] atol = atol rtol = rtol
    @test ν.basis.monomials == cert_monos

    N = SumOfSquares.Certificate.NewtonFilter{
        SumOfSquares.Certificate.NewtonDegreeBounds{Tuple{}},
    }
    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace,
        MB.SubBasis{MB.Monomial,monomial_type(x),monomial_vector_type(x)},
        SumOfSquares.Certificate.Newton{typeof(cone),MB.FullBasis{basis,monomial_type(x)},N},
    }
    @test list_of_constraint_types(model) == [(Vector{AffExpr}, S)]
    return test_delete_bridge(
        model,
        cref,
        1,
        (
            (MOI.VectorOfVariables, MOI.Nonnegatives, 0),
            (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}, 0),
            (MOI.VectorAffineFunction{Float64}, MOI.Zeros, 0),
        ),
    )
end
function sos_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), MB.Monomial, false)
end
sd_tests["sos_univariate_quadratic"] = sos_univariate_quadratic_test
function sos_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), MB.Monomial, true)
end
sd_tests["sos_bivariate_quadratic"] = sos_bivariate_quadratic_test
function sos_scaled_univariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        SOSCone(),
        MB.ScaledMonomial,
        false,
    )
end
sd_tests["sos_scaled_univariate_quadratic"] =
    sos_scaled_univariate_quadratic_test
function sos_scaled_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), MB.ScaledMonomial, true)
end
sd_tests["sos_scaled_bivariate_quadratic"] = sos_scaled_bivariate_quadratic_test
function sos_cheby_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), Chebyshev, false)
end
sd_tests["sos_cheby_univariate_quadratic"] = sos_cheby_univariate_quadratic_test
function sos_cheby_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SOSCone(), Chebyshev, true)
end
sd_tests["sos_scaled_bivariate_quadratic"] = sos_scaled_bivariate_quadratic_test

function sdsos_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), MB.Monomial, false)
end
soc_tests["sdsos_univariate_quadratic"] = sdsos_univariate_quadratic_test
function sdsos_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), MB.Monomial, true)
end
soc_tests["sdsos_bivariate_quadratic"] = sdsos_bivariate_quadratic_test
function sdsos_scaled_univariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        SDSOSCone(),
        MB.ScaledMonomial,
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
        MB.ScaledMonomial,
        true,
    )
end
soc_tests["sdsos_scaled_bivariate_quadratic"] =
    sdsos_scaled_bivariate_quadratic_test
function sdsos_cheby_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), MB.Chebyshev, false)
end
soc_tests["sdsos_cheby_univariate_quadratic"] =
    sdsos_cheby_univariate_quadratic_test
function sdsos_cheby_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, SDSOSCone(), MB.Chebyshev, true)
end
soc_tests["sdsos_scaled_bivariate_quadratic"] =
    sdsos_scaled_bivariate_quadratic_test

function dsos_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), MB.Monomial, false)
end
linear_tests["dsos_univariate_quadratic"] = dsos_univariate_quadratic_test
function dsos_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), MB.Monomial, true)
end
linear_tests["dsos_bivariate_quadratic"] = dsos_bivariate_quadratic_test
function dsos_scaled_univariate_quadratic_test(optimizer, config)
    return quadratic_test(
        optimizer,
        config,
        DSOSCone(),
        MB.ScaledMonomial,
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
        MB.ScaledMonomial,
        true,
    )
end
linear_tests["dsos_scaled_bivariate_quadratic"] =
    dsos_scaled_bivariate_quadratic_test
function dsos_cheby_univariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), MB.Chebyshev, false)
end
linear_tests["dsos_cheby_univariate_quadratic"] =
    dsos_cheby_univariate_quadratic_test
function dsos_cheby_bivariate_quadratic_test(optimizer, config)
    return quadratic_test(optimizer, config, DSOSCone(), MB.Chebyshev, true)
end
linear_tests["dsos_cheby_bivariate_quadratic"] =
    dsos_cheby_bivariate_quadratic_test
