# Adapted from:
# SOSDEMO1 --- Sum of Squares Test
# Section 3.1 of SOSTOOLS User's Manual

@testset "SOSDEMO1 with $solver" for solver in sdp_solvers
    @polyvar x y

    m = SOSModel(solver = solver)

    # Example 2.4 of
    # P. Parrilo and A. Jadbabaie
    # "Approximation of the joint spectral radius using sum of squares."
    # Linear Algebra and its Applications, Elsevier, 2008, 428, 2385-2402

    p = 2*x^4 + 2*x^3*y - x^2*y^2 + 5*y^4

    soscon = @constraint m p >= 0

    solve(m)

    @test JuMP.primalstatus(m) == MOI.FeasiblePoint

    q = getslack(soscon)
    Q = getmat(q)
    @test issymmetric(Q)
    @test isapprox(Q[1, 1], 2, rtol=1e-5)
    @test isapprox(Q[1, 2], 1, rtol=1e-5)
    @test isapprox(Q[3, 3], 5, rtol=1e-5)
    @test abs(Q[2, 3]) < 1e-5
    @test isapprox(Q[2, 2] + 2Q[1, 3], -1, rtol=1e-5)
    sosdec = SOSDecomposition(q)
    @test isapprox(sum(sosdec.ps.^2), p; rtol=1e-4, ztol=1e-6)

    M = SOSModel(solver = solver)

    p = 4*x^4*y^6 + x^2 - x*y^2 + y^2

    soscon = @constraint M p >= 0

    solve(M)

    @test JuMP.primalstatus(m) == MOI.FeasiblePoint

    # p should be MatPolynomial([1, 0, -1/2, 0, -1, 1, 0, -2/3, 0, 4/3, 0, 0, 2, 0, 4], [y, x, x*y, x*y^2, x^2*y^3])
end
