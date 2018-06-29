@testset "Motzkin with $(typeof(solver))" for solver in sdp_solvers
    @polyvar x y

    MOI.empty!(solver)
    m = SOSModel(optimizer=solver)

    p = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2

    @constraint m p >= 0

    JuMP.optimize(m)

    @test JuMP.dualstatus(m) == MOI.InfeasibilityCertificate

    MOI.empty!(solver)
    M = SOSModel(optimizer=solver)

    q = (x^2 + y^2) * p

    @constraint M q >= 0

    JuMP.optimize(M)

    @test JuMP.primalstatus(M) == MOI.FeasiblePoint
end
