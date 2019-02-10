# SPOT example from file doc/examples/sdsos_example.m
# of https://github.com/anirudhamajumdar/spotless/tree/spotless_isos

# Section 4.1 of
# A. A. Ahmadi, and A. Majumdar
# DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization
# 2017

@testset "[AM17] Section 4.1 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x[1:3]
    vx = monomials(x, 4) # Degree 4 homogeneous
    # Coefficient of polynomial
    cp = 15:-1:1
    p = polynomial(cp, vx)


    function sdsos_example(cone)
        m = SOSModel(factory)
        @variable m γ
        @constraint m p - γ*sum(x .* x)^2 in cone
        @objective m Max γ
        JuMP.optimize!(m)
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        JuMP.objective_value(m)
    end

    @test sdsos_example(DSOSCone())  ≈ -11/3     rtol=1e-5
    if !iscsdp(factory) # CSDP does not natively support SOC and uses a bridge that transforms it into SDP so it returns UNKNOWN_RESULT_STATUS on Travis
        @test sdsos_example(SDSOSCone()) ≈ -3.172412 rtol=1e-5
    end
    @test sdsos_example(SOSCone())   ≈ -0.184667 rtol=1e-5
end
