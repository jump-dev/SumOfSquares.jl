@testset "Motzkin with $(factory.constructor)" for factory in sdp_factories
    @polyvar x y

    m = SOSModel(factory)

    p = x^4*y^2 + x^2*y^4 + 1 - 3*x^2*y^2

    @constraint m p >= 0

    JuMP.optimize(m)

    @test JuMP.dualstatus(m) == MOI.InfeasibilityCertificate

    M = SOSModel(factory)

    q = (x^2 + y^2) * p

    @constraint M q >= 0

    JuMP.optimize(M)

    @test JuMP.primalstatus(M) == MOI.FeasiblePoint
end
