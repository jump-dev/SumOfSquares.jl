@testset "Polynomial equality constraint with domain with $(factory.constructor)" for factory in sdp_factories
    @polyvar(x, y)

    m = SOSModel(factory)

    @constraint(m, x^3 + x*y^2 == x, domain=(@set x^2 + y^2 >= 1 && x^2 + y^2 <= 1))

    JuMP.optimize!(m)

    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
end
