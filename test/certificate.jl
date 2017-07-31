@testset "Monomial selection for certificate" begin
    @polyvar x y
    @test getmonomialsforcertificate([x*y, y^2]) == MonomialVector([y])
    @test getmonomialsforcertificate(x^4 + 5x^2*y^2 - 2x^2*y - x*y^2 + x^2) == MonomialVector([x^2, x*y, x])
    @show getmonomialsforcertificate(x^4 + 5x^2*y^2 - 2x^2*y - x*y^2 + x^2, Polyhedra.SimplePolyhedraLibrary{Rational{BigInt}}())
    # Example 3.92
    # Blekherman, G., Parrilo, P. A., & Thomas, R. R. (Eds.).
    # Semidefinite optimization and convex algebraic geometry SIAM 2013
    @test getmonomialsforcertificate(5 - x*y - x^2*y^2 + 3y^2 + x^4) == MonomialVector([x^2, x*y, x, y, 1])
    @show getmonomialsforcertificate(5 - x*y - x^2*y^2 + 3y^2 + x^4, Polyhedra.SimplePolyhedraLibrary{Rational{BigInt}}())
end

#@testset "Random SOS should be SOS with $solver" for solver in sdp_solvers
#    @polyvar x y
#    x = [1, x, y, x^2, y^2, x*y]
#    @test_throws ArgumentError randsos(x, monotype=:Unknown)
#    for i in 1:10
#        for monotype in [:Classic, :Gram]
#            p = randsos(x, monotype=monotype)
#
#            m = SOSModel(solver = solver)
#
#            @constraint m p >= 0
#
#            status = solve(m)
#
#            @test status == :Optimal
#        end
#    end
#end
