# Adapted from:
# SOSDEMO6 --- MAX CUT
# Section 3.6 of SOSTOOLS User's Manual

@testset "SOSDEMO6 with $solver" for solver in sdp_solvers
    @polyvar x[1:5]

    # Number of cuts
    f = 2.5 - 0.5*x[1]*x[2] - 0.5*x[2]*x[3] - 0.5*x[3]*x[4] - 0.5*x[4]*x[5] - 0.5*x[5]*x[1]

    # Boolean constraints
    bc = vec(x).^2 - 1

    for (gamma, expected) in [(3.9, MOI.InfeasiblePoint), (4, MOI.FeasiblePoint)]

        m = SOSModel(solver = solver)

        Z = monomials(x, 0:1)
        @variable m p1 SOSPoly(Z)

        Z = monomials(x, 0:2)
        @variable m p[1:5] Poly(Z)

        @constraint m p1*(gamma-f) + dot(p, bc) >= (gamma-f)^2

        solve(m)

        @test JuMP.primalstatus(m) == expected
    end
end
