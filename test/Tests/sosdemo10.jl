# Adapted from:
# SOSDEMO10 --- Set containment
# Section 3.10 of SOSTOOLS User's Manual

@testset "SOSDEMO10 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x[1:2]

    eps = 1e-6
    p = x[1]^2+x[2]^2
    gamma = 1
    g0 = 2*x[1]
    theta = 1

    m = SOSModel(factory)

    # FIXME s should be sos ?
    # in SOSTools doc it is said to be SOS
    # but in the demo it is not constrained so
    Z = monomials(x, 0:4)
    @variable m s Poly(Z)

    Z = monomials(x, 2:3)
    @variable m g1 Poly(Z)

    Sc = [theta^2-s*(gamma-p) g0+g1; g0+g1 1]

    @SDconstraint m Matrix(eps * I, 2, 2) ⪯ Sc

    JuMP.optimize!(m)

    # Program is feasible, that is, the set
    # { x |((g0+g1) + theta)(theta - (g0+g1)) >=0 }
    # contains the set
    # { x | p <= gamma }

    # ALMOST_OPTIMAL and NEARLY_FEASIBLE_POINT for CSDP
    @test JuMP.termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
    @test JuMP.primal_status(m) in [MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT]

    #@show JuMP.value(s)
    #@show JuMP.value(g1)
end
