# Example 3.25, 3.35 and 3.38 of
# Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
# Semidefinite optimization and convex algebraic geometry SIAM 2013
@testset "Example 3.25 with $solver" for solver in sdp_solvers
    @polyvar x y
    m = Model(solver=solver)
    @polyconstraint m x^2-x*y^2+y^4+1 >= 0
    status = solve(m)
    @test status == :Optimal
end
@testset "Example 3.35 with $solver" for solver in sdp_solvers
    @polyvar x
    m = Model(solver=solver)
    @polyconstraint m x^4+4x^3+6x^2+4x+5 >= 0
    status = solve(m)
    @test status == :Optimal
    # p(x) = (x^2+2x)^2 + 2(1+x)^2 + 3
end
@testset "Example 3.38 with $solver" for solver in sdp_solvers
    @polyvar x y
    m = Model(solver=solver)
    @polyconstraint m 2x^4 + 5y^4 - x^2*y^2 >= -2(x^3*y + x + 1)
    status = solve(m)
    @test status == :Optimal
    # p(x) = (x^2+2x)^2 + 2(1+x)^2 + 3
end
