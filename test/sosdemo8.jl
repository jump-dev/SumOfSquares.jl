# Adapted from:
# SOSDEMO8 --- Bounds in Probability
# Section 3.8 of SOSTOOLS User's Manual

@testset "SOSDEMO8 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x

    # The probability adds up to one.
    μ0 = 1
    # Mean
    μ1  = 1
    # Variance
    σ = 1/2
    # E(x^2)
    μ2 = σ^2+μ1^2
    # Support of the random variable
    R = [0,5]

    # Event whose probability we want to bound
    E = [4,5]

    m = SOSModel(factory)

    @variable m a
    @variable m b
    @variable m c

    P = a + b*x + c*x^2

    # Nonnegative on the support
    @constraint(m, P >= 0, domain = (@set 0 <= x && x <= 5))

    # Greater than one on the event
    @constraint(m, P >= 1, domain = (@set 4 <= x && x <= 5))

    # The bound
    bnd =  a * μ0 + b * μ1 + c * μ2

    @objective m Min bnd

    JuMP.optimize!(m)
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    @test isapprox(JuMP.objective_value(m), 1/37, rtol=1e-5)
    @test isapprox(JuMP.value(P), ((12/37)x-11/37)^2, rtol=1e-3)
end
