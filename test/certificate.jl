@testset "Monomial selection for certificate" begin
    @polyvar x y z
    @testset "Multipartite error not disjoint: $parts" for parts in [([x], [x]),
                                                                     ([x], [x, y]),
                                                                     ([x], [y], [x, y])]
        err = ArgumentError("Parts are not disjoint in multipartite Newton polytope estimation: $parts")
        @test_throws err SumOfSquares.monomials_half_newton_polytope([x*y, y^2], parts)
    end
    @testset "Unipartite" begin
        @test SumOfSquares.monomials_half_newton_polytope([x*y, y^2], tuple()) == [y]
        @test isempty(SumOfSquares.monomials_half_newton_polytope([x, y], tuple()))
        @test SumOfSquares.monomials_half_newton_polytope([x^2, y^2], tuple()) == [x, y]
        @test SumOfSquares.monomials_half_newton_polytope([x^2, y^2], ([x, y],)) == [x, y]
        @test SumOfSquares.monomials_half_newton_polytope([x^2, y^2], ([y, x],)) == [x, y]
    end
    @testset "Multipartite" begin
        # In the part [y, z], the degree is between 0 and 2
        X = [x^4, x^2*y^2, x^2*z^2, x^2*y*z, y*z]
        @test SumOfSquares.monomials_half_newton_polytope(X, tuple()) == [x^2, x*y, x*z, y*z, x, y, z]
        function full_test(X, Y, part1, part2)
            @test SumOfSquares.monomials_half_newton_polytope(X, (part1,)) == Y
            @test SumOfSquares.monomials_half_newton_polytope(X, (part2,)) == Y
            a = SumOfSquares.monomials_half_newton_polytope(X, (part2,))
            @test SumOfSquares.monomials_half_newton_polytope(X, (part1, part2)) == Y
            @test SumOfSquares.monomials_half_newton_polytope(X, (part2, part1)) == Y
        end
        full_test(X, monovec([x^2, x*y, x*z,      x, y, z]), [x], [y, z])
        full_test(X, monovec([x^2, x*y, x*z, y*z, x,    z]), [y], [x, z])
        full_test(X, monovec([x^2, x*y, x*z, y*z, x, y   ]), [z], [x, y])
        # FIXME: With recursive merging, it should give [x^2, x*y, x*z, x]
        @test SumOfSquares.monomials_half_newton_polytope([x^4, x^2*y^2, x^2*z^2, x^2*y*z, y*z], ([x], [y], [z])) == [x^2, x*y, x*z, y*z, x, y, z]
    end
end

@testset "Random SOS should be SOS with $(factory.constructor)" for factory in sdp_factories
    @polyvar x y
    x = [1, x, y, x^2, y^2, x*y]
    @test_throws ArgumentError randsos(x, monotype=:Unknown)
    for i in 1:10
        for monotype in [:Classic, :Gram]
            p = randsos(x, monotype=monotype)

            m = SOSModel(factory)

            @constraint m p >= 0

            JuMP.optimize!(m)

            @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        end
    end
end
