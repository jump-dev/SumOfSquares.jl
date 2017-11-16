# Adapted from:
# SOSDEMO9 --- Matrix SOS decomposition
# Section 3.9 of SOSTOOLS User's Manual

@testset "SOSDEMO9 with $solver" for solver in sdp_solvers
    @polyvar x1 x2 x3

    P = [x1^4+x1^2*x2^2+x1^2*x3^2 x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2);
         x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2) x1^2*x2^2+x2^2*x3^2+(x2^2+2*x3^2)^2]

    # Test if P(x1,x2,x3) is an SOS matrix
    m = SOSModel(solver = solver)
    # TODO return H so that P = H.'*H
    @SDconstraint m P ⪰ 0

    solve(m)
    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
end
