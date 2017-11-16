# Adapted from:
# SOSDEMO10 --- Set containment
# Section 3.10 of SOSTOOLS User's Manual

@testset "SOSDEMO10 with $solver" for solver in sdp_solvers
    @polyvar x[1:2]

    eps = 1e-6
    p = x[1]^2+x[2]^2
    gamma = 1
    g0 = 2*x[1]
    theta = 1

    m = SOSModel(solver = solver)

    # FIXME s should be sos ?
    # in SOSTools doc it is said to be SOS
    # but in the demo it is not constrained so
    Z = monomials(x, 0:4)
    @variable m s Poly(Z)

    Z = monomials(x, 2:3)
    @variable m g1 Poly(Z)

    Sc = [theta^2-s*(gamma-p) g0+g1; g0+g1 1]

    @SDconstraint m eps * eye(2) ⪯ Sc

    solve(m)

    # Program is feasible, { x |((g0+g1) + theta)(theta - (g0+g1)) >=0 } contains { x | p <= gamma }
    @test JuMP.primalstatus(m) == MOI.FeasiblePoint || JuMP.primalstatus(m) == MOI.NearlyFeasiblePoint

    #@show JuMP.resultvalue(s)
    #@show JuMP.resultvalue(g1)
end
