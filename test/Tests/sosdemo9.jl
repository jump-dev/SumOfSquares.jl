# Adapted from:
# SOSDEMO9 --- Matrix SOS decomposition
# Section 3.9 of SOSTOOLS User's Manual

@testset "SOSDEMO9 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x1 x2 x3

    P = [x1^4+x1^2*x2^2+x1^2*x3^2 x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2);
         x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2) x1^2*x2^2+x2^2*x3^2+(x2^2+2*x3^2)^2]

    # Test if P(x1,x2,x3) is an SOS matrix
    m = SOSModel(factory)
    # TODO return H so that P = H.'*H
    @SDconstraint m P ⪰ 0

    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    # NEARLY_FEASIBLE_POINT for CSDP
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
end
