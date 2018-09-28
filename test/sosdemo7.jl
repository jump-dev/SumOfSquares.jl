# Adapted from:
# SOSDEMO7 --- Chebyshev polynomials
# Section 3.7 of SOSTOOLS User's Manual

@testset "SOSDEMO7 with $(factory.constructor)" for factory in sdp_factories
    isscs(factory) && continue
    ndeg = 8   # Degree of Chebyshev polynomial

    @polyvar x

    Z = monomials((x,), 0:ndeg-1)

    m = SOSModel(factory)

    @variable m γ
    @variable m p1 Poly(Z)

    p = p1 + γ * x^ndeg # the leading coeff of p is γ

    dom = @set x >= -1 && x <= 1
    @constraint(m, p <= 1, domain = dom)
    @constraint(m, p >= -1, domain = dom)

    @objective m Max γ

    JuMP.optimize!(m)

    @test JuMP.primal_status(m) == MOI.FeasiblePoint

    @test isapprox(JuMP.result_value(p), 128x^8 - 256x^6 + 160x^4 - 32x^2 + 1, ztol=1e-7, atol=1e-7)
    @test isapprox(JuMP.result_value(γ), 128)
end
