# SPOT example from file doc/examples/sdsos_example.m
# of https://github.com/anirudhamajumdar/spotless/tree/spotless_isos

# Section 4.1 of
# A. A. Ahmadi, and A. Majumdar
# DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization
# 2017

@testset "[AM17] Section 4.1 with $solver" for solver in sdp_solvers
    @polyvar x[1:3]
    vx = monomials(x, 4) # Degree 4 homogeneous
    # Coefficient of polynomial
    cp = collect(15:-1:1) # TODO remove collect once DynamicPolynomials is released
    p = polynomial(cp, vx)


    function sdsos_example(cone)
        m = SOSModel(solver = solver)
        @variable m γ
        @constraint m p - γ*sum(x .* x)^2 in cone
        @objective m Max γ
        status = solve(m)
        @test status == :Optimal
        getobjectivevalue(m)
    end

    @test sdsos_example(DSOSCone()) ≈ -11/3 rtol=1e-5
    @test sdsos_example(SOSCone()) ≈ -0.184667 rtol=1e-5
end
