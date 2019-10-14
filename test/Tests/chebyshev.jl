# Adapted from:
# SOSDEMO7 --- Chebyshev polynomials
# Section 3.7 of SOSTOOLS User's Manual

using Test
using SumOfSquares
using DynamicPolynomials

function chebyshev_test(optimizer, config::MOIT.TestConfig)
    ndeg = 8   # Degree of Chebyshev polynomial

    @polyvar x

    Z = monomials((x,), 0:ndeg-1)

    model = _model(optimizer)

    @variable(model, γ)
    @variable(model, p1, Poly(Z))

    p = p1 + γ * x^ndeg # the leading coeff of p is γ

    dom = @set x >= -1 && x <= 1
    @constraint(model, 1 - p in SOSCone(), domain = dom)
    @constraint(model, p + 1 in SOSCone(), domain = dom)

    @objective(model, Max, γ)

    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT

    @test isapprox(JuMP.value(p), 128x^8 - 256x^6 + 160x^4 - 32x^2 + 1, ztol=config.atol, atol=config.atol, rtol=config.rtol)
    @test isapprox(JuMP.value(γ), 128, atol=config.atol, rtol=config.rtol)
end
sd_tests["chebyshev"] = chebyshev_test
