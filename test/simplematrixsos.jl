# Example 3.77 and 3.79 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
@testset "Example 3.77 and 3.79 with $(typeof(solver))" for solver in sdp_solvers
    @polyvar x
    P = [x^2-2x+2 x; x x^2]
    # Example 3.77
    MOI.empty!(solver)
    m = SOSModel(optimizer=solver)
    @SDconstraint m P >= 0
    JuMP.optimize(m)
    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
    # Example 3.79
    @polyvar y[1:2]
    MOI.empty!(solver)
    M = SOSModel(optimizer=solver)
    @constraint M dot(vec(y), P*vec(y)) >= 0
    JuMP.optimize(M)
    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
end
