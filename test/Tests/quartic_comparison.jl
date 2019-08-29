# SPOT example from file doc/examples/sdsos_example.m
# of https://github.com/anirudhamajumdar/spotless/tree/spotless_isos
# See Section 4.1 of [AM17].
#
# [AM17] Section 4.1 of
# A. A. Ahmadi, and A. Majumdar
# DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization
# 2017

using Test
using SumOfSquares
using DynamicPolynomials

function quartic_comparison_test(
    optimizer, config::MOIT.TestConfig,
    cone::SumOfSquares.PolyJuMP.PolynomialSet, expected_objective_value)
    atol = config.atol
    rtol = config.rtol

    @polyvar x[1:3]
    vx = monomials(x, 4) # Degree 4 homogeneous
    # Coefficient of polynomial
    cp = 15:-1:1
    p = polynomial(cp, vx)

    model = _model(optimizer)
    @variable(model, γ)
    @constraint(model, p - γ * sum(x .* x)^2 in cone)
    @objective(model, Max, γ)
    JuMP.optimize!(model)
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    @test JuMP.objective_value(model) ≈ expected_objective_value atol=atol rtol=rtol
end

sos_quartic_comparison_test(optimizer, config)   = quartic_comparison_test(optimizer, config, SOSCone(), -0.184667)
sd_tests["sos_quartic_comparison"]      = sos_quartic_comparison_test
sdsos_quartic_comparison_test(optimizer, config) = quartic_comparison_test(optimizer, config, SDSOSCone(), -3.172412)
soc_tests["sdsos_quartic_comparison"]   = sdsos_quartic_comparison_test
dsos_quartic_comparison_test(optimizer, config)  = quartic_comparison_test(optimizer, config, DSOSCone(), -11/3)
linear_tests["dsos_quartic_comparison"] = dsos_quartic_comparison_test
