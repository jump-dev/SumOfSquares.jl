# Example 3.77 and 3.79 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
@testset "Example 3.77 and 3.79 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x
    P = [x^2-2x+2 x; x x^2]
    # Example 3.77
    m = SOSModel(factory)
    @SDconstraint m P >= 0
    JuMP.optimize(m)
    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
    # Example 3.79
    @polyvar y[1:2]
    M = SOSModel(factory)
    @constraint M dot(vec(y), P*vec(y)) >= 0
    JuMP.optimize(M)
    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
end
