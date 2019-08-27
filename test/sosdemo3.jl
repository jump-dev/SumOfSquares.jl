# Adapted from:
# SOSDEMO3 --- Bound on Global Extremum
# Section 3.3 of SOSTOOLS User's Manual

@testset "SOSDEMO3 with $(factory.constructor)" for factory in sdp_factories
    iscsdp(factory) && continue # will be fixed by objective bridge
    @polyvar x1 x2

    m = SOSModel(factory)

    @variable m γ

    # Constraint : r(x)*(f(x) - gam) >= 0
    # f(x) is the Goldstein-Price function
    f1 = x1+x2+1
    f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
    f3 = 2*x1-3*x2
    f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2

    f = (1+f1^2*f2)*(30+f3^2*f4)

    @constraint m f >= γ

    @objective m Max γ

    JuMP.optimize!(m)

    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT || JuMP.primal_status(m) == MOI.NEARLY_FEASIBLE_POINT

    @test isapprox(JuMP.objective_value(m), 3; rtol=1e-2)
end
