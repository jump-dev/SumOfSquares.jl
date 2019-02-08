# Example 3.77 and 3.79 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
@testset "Example 3.77 and 3.79 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x
    P = [x^2-2x+2 x; x x^2]
    # Example 3.77
    m = SOSModel(factory)
    mat_cref = @SDconstraint m P >= 0
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test length(gram_matrix(mat_cref).x) == 4
    # Example 3.79
    @polyvar y[1:2]
    M = SOSModel(factory)
    cref = @constraint M dot(vec(y), P*vec(y)) >= 0
    JuMP.optimize!(M)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test length(gram_matrix(cref).x) == 6
end
