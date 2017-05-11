@testset "Polynomial equality constraint with domain" for solver in sdp_solvers
    @polyvar(x, y)

    m = SOSModel(solver=solver)

    @polyconstraint(m, x^3 + x*y^2 == x, domain=(x^2 + y^2 >= 1 && x^2 + y^2 <= 1))

    # @variable(m, γ)
    # @polyconstraint(m, γ >= x + y)

    # @objective(m, Min, γ)

    status = solve(m)

    @test status == :Optimal || (iscsdp(solver) && status == :Suboptimal)

    @test isapprox(getobjectivevalue(m), 3; rtol=1e-3)
end
