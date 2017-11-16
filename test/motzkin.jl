@testset "Motzkin with $solver" for solver in sdp_solvers
    @polyvar x y

    m = SOSModel(solver = solver)

    p = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2

    @constraint m p >= 0

    solve(m)

    @test JuMP.primalstatus(m) == MOI.InfeasiblePoint

    M = SOSModel(solver = solver)

    q = (x^2 + y^2) * p

    @constraint M q >= 0

    solve(M)

    @test JuMP.primalstatus(M) == MOI.FeasiblePoint
end
