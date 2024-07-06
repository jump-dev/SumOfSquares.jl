# Section 4.4 of
# A. A. Ahmadi, and A. Majumdar
# DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization
# 2017

using MultivariateMoments

function options_pricing_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
    K::Int,
    expected::Float64,
)
    VERSION < v"1.0.3" && return # see https://github.com/jump-dev/SumOfSquares.jl/issues/48
    atol = 100config.atol
    rtol = 100config.rtol

    @polyvar x y z
    σ = [184.04, 164.88, 164.88, 184.04, 164.88, 184.04]
    X = [x^2, x * y, x * z, y^2, y * z, z^2, x, y, z, 1]
    μ = moment_vector([σ .+ 44.21^2; 44.21 * ones(3); 1], X)

    cocone = CopositiveInner(cone)

    model = _model(optimizer)
    @variable(model, p, Poly(X))
    @constraint(model, p in cocone)
    @constraint(model, p - (x - K) in cocone)
    @constraint(model, p - (y - K) in cocone)
    @constraint(model, p - (z - K) in cocone)
    @objective(model, Min, dot(μ, p))
    if MOI.Test._supports(config, MOI.optimize!)
        JuMP.optimize!(model)
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        @test JuMP.objective_value(model) ≈ expected atol = atol rtol = rtol
    end
end

function options_pricing_test(
    optimizer,
    config::MOI.Test.Config,
    cone::SumOfSquares.PolyJuMP.PolynomialSet,
    Ks::Vector{Int},
    expected::Vector{Float64},
)
    @testset "K = $K" for (K, exp) in zip(Ks, expected)
        options_pricing_test(optimizer, config, cone, K, exp)
    end
end

const K = [30, 35, 40, 45, 50]
const dsos_codsos_exp = [132.63, 132.63, 132.63, 132.63, 132.63]
const sdsos_cosdsos_exp = [21.51, 17.17, 13.20, 9.85, 7.30]

function sos_options_pricing_test(optimizer, config)
    return options_pricing_test(
        optimizer,
        config,
        SOSCone(),
        K,
        sdsos_cosdsos_exp,
    )
end
sd_tests["sos_options_pricing"] = sos_options_pricing_test
function sdsos_options_pricing_test(optimizer, config)
    return options_pricing_test(
        optimizer,
        config,
        SDSOSCone(),
        K,
        sdsos_cosdsos_exp,
    )
end
soc_tests["sdsos_options_pricing"] = sdsos_options_pricing_test
function dsos_options_pricing_test(optimizer, config)
    return options_pricing_test(
        optimizer,
        config,
        DSOSCone(),
        K,
        dsos_codsos_exp,
    )
end
linear_tests["dsos_options_pricing"] = dsos_options_pricing_test
