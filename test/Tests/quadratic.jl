using Test
import MultivariateBases
using DynamicPolynomials

function _test_moments(test_values, μ, monos)
    @test μ isa AbstractMeasure{Float64}
    @test length(moments(μ)) == length(monos)
    test_values(moment_value.(moments(μ)))
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
        if basis === Chebyshev
            # See https://github.com/jump-dev/SumOfSquares.jl/issues/357
            cert_monos = [1, y, x]
        else
            cert_monos = [y, x]
        end
        if basis === Chebyshev
            monos = [1, y^2, x * y, x^2]
        else
            monos = [y^2, x * y, x^2]
        end
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

    test_constraint_primal(cref, value(poly); atol, rtol)

    p = gram_matrix(cref)
    if basis === Chebyshev && bivariate
        # See https://github.com/jump-dev/SumOfSquares.jl/issues/357
        @test value_matrix(p) ≈ [zeros(1, 3); zeros(2) ones(2, 2)] atol = atol rtol =
            rtol
    else
        @test value_matrix(p) ≈ ones(2, 2) atol = atol rtol = rtol
    end
    @test p.basis.monomials == cert_monos

    μ = moments(dual(cref))
    a = moment_value.(μ)
    if bivariate && basis === MB.Chebyshev
        b = a[2:end] + [a[1], 0, a[1]]
        b[1] /= 2
        b[3] /= 2
    else
        b = a
    end
    @test b[2] ≈ (bivariate && basis === MB.ScaledMonomial ? -√2 : -1.0) atol =
        atol rtol = rtol
    @test b[1] + b[3] ≈ 2.0 atol = atol rtol = rtol
    @test μ[2].polynomial == MB.Polynomial{basis}(
        bivariate ? (basis === MB.Chebyshev ? y^2 : x * y) : x^1,
    )

    @test dual_status(model) == MOI.FEASIBLE_POINT
    _test_moments(dual(cref), monos) do vals
        @test vals ≈ a atol = atol rtol = rtol
    end
    if basis === MB.Chebyshev && bivariate
        _test_moments(
            moments(cref),
            monomial_vector([1, y, x, y^2, x * y, x^2]),
        ) do vals
            @test length(vals) == 6
            for i in [2, 3]
                @test vals[i] ≈ 0 rtol = rtol atol = atol
            end
            @test vals[5] ≈ -1 rtol = rtol atol = atol
            @test vals[4] + vals[6] + 2vals[1] ≈ 4 rtol = rtol atol = atol
        end
    else
        _test_moments(moments(cref), monos) do vals
            @test vals ≈ a atol = atol rtol = rtol
        end
    end

    ν = moment_matrix(cref)
    off = if bivariate && basis === ScaledMonomial
        a[2] / √2
    elseif bivariate && basis == Chebyshev
        -(a[1] + a[2]) / 2
    else
        a[2]
    end
    M = value_matrix(ν)
    if basis === Chebyshev && bivariate
        M = M[2:end, 2:end]
    end
    @test M ≈ [
        b[1] off
        off b[3]
    ] atol = atol rtol = rtol
    @test ν.basis.monomials == cert_monos

    N = SumOfSquares.Certificate.NewtonFilter{
        SumOfSquares.Certificate.NewtonDegreeBounds{Tuple{}},
    }
    S = SumOfSquares.SOSPolynomialSet{
        SumOfSquares.FullSpace,
        MB.SubBasis{basis,monomial_type(x),monomial_vector_type(x)},
        SumOfSquares.Certificate.Newton{
            typeof(cone),
            MB.FullBasis{basis,monomial_type(x)},
            MB.FullBasis{basis,monomial_type(x)},
            N,
        },
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
