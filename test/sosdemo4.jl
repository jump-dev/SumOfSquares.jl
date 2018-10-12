# Adapted from:
# SOSDEMO4 --- Matrix Copositivity
# Section 3.4 of SOSTOOLS User's Manual

@testset "SOSDEMO4 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x[1:5]

    # The matrix under consideration
    J = [1 -1  1  1 -1;
        -1  1 -1  1  1;
         1 -1  1 -1  1;
         1  1 -1  1 -1;
        -1  1  1 -1  1]

    xs = vec(x).^2
    xsJxs = dot(xs, J*xs)
    r = sum(xs)

    m0 = SOSModel(factory)
    @constraint m0 xsJxs >= 0
    JuMP.optimize!(m0)
    @test JuMP.dual_status(m0) == MOI.InfeasibilityCertificate

    m1 = SOSModel(factory)
    @constraint m1 r*xsJxs >= 0
    JuMP.optimize!(m1)
    @test JuMP.primal_status(m1) == MOI.FeasiblePoint
    # Program is feasible. The matrix J is copositive
end
