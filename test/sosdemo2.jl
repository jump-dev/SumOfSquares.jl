# Adapted from:
# SOSDEMO2 --- Lyapunov Function Search
# Section 3.2 of SOSTOOLS User's Manual

@testset "SOSDEMO2 with $solver" for solver in sdp_solvers
    @polyvar x[1:3]

    # Constructing the vector field dx/dt = f
    f = [-x[1]^3-x[1]*x[3]^2,
         -x[2]-x[1]^2*x[2],
    -x[3]+3*x[1]^2*x[3]-3*x[3]/(x[3]^2+1)]

    m = SOSModel(solver = solver)

    # The Lyapunov function V(x):
    Z = vec(x).^2
    @variable m V Poly(Z)

    @constraint m V >= sum(x.^2)

    # dV/dx*(x[3]^2+1)*f <= 0
    P = dot(differentiate(V, x), f).num # the denominator is x[3]^2+1
    @constraint m P <= 0

    solve(m)

    @test JuMP.primalstatus(m) == MOI.FeasiblePoint

    @test iszero(removemonomials(JuMP.resultvalue(V), Z))
end
