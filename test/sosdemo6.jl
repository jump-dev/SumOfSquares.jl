# Adapted from:
# SOSDEMO6 --- MAX CUT
# Section 3.6 of SOSTOOLS User's Manual

@testset "SOSDEMO6 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x[1:5]

    # Number of cuts
    f = 2.5 - 0.5*x[1]*x[2] - 0.5*x[2]*x[3] - 0.5*x[3]*x[4] - 0.5*x[4]*x[5] - 0.5*x[5]*x[1]

    # Boolean constraints
    bc = vec(x).^2 - 1

    for (gamma, feasible) in [(3.9, false), (4, true)]

        m = SOSModel(factory)

        Z = monomials(x, 0:1)
        @variable m p1 SOSPoly(Z)

        Z = monomials(x, 0:2)
        @variable m p[1:5] Poly(Z)

        @constraint m p1*(gamma-f) + dot(p, bc) >= (gamma-f)^2

        JuMP.optimize!(m)

        if feasible
            @test JuMP.primal_status(m) == MOI.FeasiblePoint
        else
            @test JuMP.dual_status(m) == MOI.InfeasibilityCertificate
        end
    end
end
