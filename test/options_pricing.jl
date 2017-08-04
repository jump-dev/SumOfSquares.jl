# Section 4.4 of
# A. A. Ahmadi, and A. Majumdar
# DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization
# 2017

using MultivariateMoments

@testset "Options Pricing with $solver" for solver in sdp_solvers
    isscs(solver) && continue
    @polyvar x y z
    σ = [184.04, 164.88, 164.88, 184.04, 164.88, 184.04]
    X2 = [x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1]
    X = [x, y, z, 1]
    μ = Measure([σ + 44.21^2; 44.21 * ones(3); 1],
                X2)
    function optionspricing(K, poly, cone)
        m = SOSModel(solver = solver)
        @variable m p poly(X)
        @constraint m p - (x - K) in cone
        @constraint m p - (y - K) in cone
        @constraint m p - (z - K) in cone
        @objective m Min dot(μ, p)
        status = solve(m)
        @test status == :Optimal
        getobjectivevalue(m)
    end

    K = [30, 35, 40, 45, 50]
    #dsos_codsos_exp = [132.63, 132.63, 132.63, 132.63, 132.63]
    dsos_codsos_exp = [1670.99, 1670.99, 1670.99, 1670.99, 1670.99]
    dsos_cosos_exp = [21.73, 18.53, 16.16, 14.32, 12.86]
    sos_codsos_exp = [265.50, 265.50, 265.50, 265.50, 265.50]
    sos_cosos_exp = [21.51, 17.17, 13.20, 9.85, 7.30]
    for i in eachindex(K)
        @test optionspricing(K[i], DSOSPoly, CoDSOSCone()) ≈ dsos_codsos_exp[i] atol=1e-2
        @test optionspricing(K[i], DSOSPoly, CoSOSCone()) ≈ dsos_cosos_exp[i] atol=1e-2
        @test optionspricing(K[i], SOSPoly, CoDSOSCone()) ≈ sos_codsos_exp[i] atol=1e-2
        @test optionspricing(K[i], SOSPoly, CoSOSCone()) ≈ sos_cosos_exp[i] atol=1e-2
    end
end
