# Adapted from:
# SOSDEMO1 --- Sum of Squares Test
# Section 3.1 of SOSTOOLS User's Manual

@testset "SOSDEMO1 with $(factory.constructor)" for factory in sdp_factories
    @polyvar x y

    m = SOSModel(factory)

    # Example 2.4 of
    # P. Parrilo and A. Jadbabaie
    # "Approximation of the joint spectral radius using sum of squares."
    # Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402

    p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4

    soscon = @constraint m p >= 0

    JuMP.optimize!(m)

    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    q = gram_matrix(soscon)
    Q = getmat(q)
    @test issymmetric(Q)
    @test isapprox(Q[1, 1], 2, rtol=1e-5)
    @test isapprox(Q[1, 2], 1, rtol=1e-5)
    @test isapprox(Q[3, 3], 5, rtol=1e-5)
    @test abs(Q[2, 3]) < 1e-5
    @test isapprox(Q[2, 2] + 2Q[1, 3], -1, rtol=1e-5)
    sosdec = SOSDecomposition(q)
    @test isapprox(sosdec, sos_decomposition(soscon))
    @test isapprox(sum(sosdec.ps.^2), p; rtol=1e-4, ztol=1e-6)

    M = SOSModel(factory)

    p = 4*x^4*y^6 + x^2 - x*y^2 + y^2

    soscon = @constraint M p >= 0

    JuMP.optimize!(M)

    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT

    # p should be GramMatrix([1, 0, -1/2, 0, -1, 1, 0, -2/3, 0, 4/3, 0, 0, 2, 0, 4], [y, x, x*y, x*y^2, x^2*y^3])
end
