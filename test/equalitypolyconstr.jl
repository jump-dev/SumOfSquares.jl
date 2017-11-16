@testset "Polynomial equality constraint with domain" for solver in sdp_solvers
    @polyvar(x, y)

    m = SOSModel(solver=solver)

    @constraint(m, x^3 + x*y^2 == x, domain=(@set x^2 + y^2 >= 1 && x^2 + y^2 <= 1))

    solve(m)

    @test JuMP.primalstatus(m) == MOI.FeasiblePoint
end
