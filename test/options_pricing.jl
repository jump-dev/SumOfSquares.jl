# Section 4.4 of
# A. A. Ahmadi, and A. Majumdar
# DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization
# 2017

using MultivariateMoments

@testset "Options Pricing with $(typeof(solver))" for solver in sdp_solvers
    isscs(solver) && continue
    @polyvar x y z
    σ = [184.04, 164.88, 164.88, 184.04, 164.88, 184.04]
    X = [x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1]
    μ = measure([σ + 44.21^2; 44.21 * ones(3); 1],
                X)
    function optionspricing(K, cone)
        MOI.empty!(solver)
        m = SOSModel(optimizer=solver)
        @variable m p Poly(X)
        @constraint m p in cone
        @constraint m p - (x - K) in cone
        @constraint m p - (y - K) in cone
        @constraint m p - (z - K) in cone
        @objective m Min dot(μ, p)
        JuMP.optimize(m)
        @test JuMP.primalstatus(m) == MOI.FeasiblePoint
        JuMP.objectivevalue(m)
    end

    K = [30, 35, 40, 45, 50]
    dsos_codsos_exp = [132.63, 132.63, 132.63, 132.63, 132.63]
    sos_cosos_exp = [21.51, 17.17, 13.20, 9.85, 7.30]
    for i in eachindex(K)
        @test optionspricing(K[i], CoDSOSCone()) ≈ dsos_codsos_exp[i] atol=1e-2
        @test optionspricing(K[i], CoSOSCone()) ≈ sos_cosos_exp[i] atol=1e-2
    end
end
