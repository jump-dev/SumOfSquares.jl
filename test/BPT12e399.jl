# Adapted from:
# Example 3.99 of
# Blekherman, G.; Parrilo, P. & Thomas, R.
# Semidefinite Optimization and Convex Algebraic Geometry
# Society for Industrial and Applied Mathematics, 2012

@testset "[BPT12] Example 3.99 with $solver" for solver in sdp_solvers
    @polyvar x y

    m = SOSModel(solver = solver)

    @variable m α

    @constraint(m, x^2 + α*y <= 10, domain = @set x^2 + y^2 == 1)

    @objective m Max α

    solve(m)

    @test JuMP.primalstatus(m) == MOI.FeasiblePoint

    @test JuMP.resultvalue(α) ≈ 6 atol=1e-6

    # FIXME JuMP does not support modification once loaded on jump/moi yet
#    @objective m Min α
#
#    solve(m)
#
#    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
#
#    @test JuMP.resultvalue(α) ≈ -6 atol=1e-6
end
