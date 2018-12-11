@testset "MatPolynomial tests" begin
    @testset "MatPolynomial" begin
        @polyvar x y
        P = MatPolynomial{Int}((i,j) -> i + j, [x^2, x*y, y^2])
        @test variables(P) == [x, y]
        zP = zero(typeof(P))
        @test isempty(zP.Q)
        @test zP == 0
        p = polynomial(P)
        @test coefficients(p) == [2, 6, 12, 10, 6]
        @test monomials(p) == [x^4, x^3*y, x^2*y^2, x*y^3, y^4]
        for i in 1:3
            for j in 1:3
                @test P[i, j] == i + j
            end
        end
        for P in (MatPolynomial((i,j) -> i * j, [y, x]),
                  MatPolynomial((i,j) -> (3-i) * (3-j), monovec([y, x])),
                  MatPolynomial([1 2; 2 4], [y, x]),
                  MatPolynomial([4 2; 2 1], monovec([y, x])))
            @test P.Q.Q == [4, 2, 1]
            @test P.x[1] == x
            @test P.x[2] == y
        end
        P = MatPolynomial((i,j) -> ((i,j) == (1,1) ? 2 : 0), [x*y, x^2, y^2])
        Q = MatPolynomial([0 1; 1 0], [x^2, y^2])
        @test P == Q
        p = MatPolynomial([2 3; 3 2], [x, y])
        @test polynomial(p) isa AbstractPolynomial{Int}
        @test polynomial(p, Float64) isa AbstractPolynomial{Float64}

        @test -2*x + dot(-x - x^2, 0) + MatPolynomial{Int}((i,j)->1, [1,x]) == -(-x^2 - 1)
        P = MatPolynomial{Int}((i,j) -> i + j, [x^2, x*y, y^2])
        p = polynomial(P)
        @test !iszero(P)
        @test iszero(P-P)
        @test iszero(P-p)
        @test iszero(p-P)
        @test P + P == P + p
        @test x * P == x * p
        @test 2 * P == 2 * p
        @test P * (2x) == (2x) * p

        @test differentiate(MatPolynomial{Int}((i,j)->1, [x]), y) == 0
        @test differentiate(MatPolynomial{Int}((i,j)->1, [x, y]), y, 1) == 2x + 2y
        #@inferred differentiate(MatPolynomial{Int}((i,j)->1, [x, y]), y, 0) # FIXME failing at the moment

        @test MatPolynomial([2 3; 3 2], [x, y]) == 2x^2 + 2y^2 + 6x*y

        mp = MatPolynomial{Int}((i,j) -> i+j, [x])
        @test norm(mp) == norm(mp, 2) == 2.0

        @polyvar v[1:3]
        P = [1 2 3; 2 4 5; 3 5 6]
        p = MatPolynomial(P, v)
        #@inferred p(x => ones(3))
        @test p(v => ones(3)) == 31
        #@inferred subs(p, x => ones(3))
        @test subs(p, v => ones(3)) == 31
    end
    @testset "SOSDecomposition" begin
        @polyvar x y
#       @test isempty(SOSDecomposition(typeof(x)[]))
#       ps = [1, x + y, x^2, x*y, 1 + x + x^2]
#       P = MatPolynomial(SOSDecomposition(ps))
#       P.Q == [2 0 1 0 1; 0 1 0 0 0; 1 0 2 1 1; 0 0 1 1 0; 1 0 1 0 2]
#       P.x == [x^2, x*y, x, y, 1]
#       @test P == P
#       @test isapprox(MatPolynomial(SOSDecomposition(P)), P)
        @test sprint(show, SOSDecomposition([x+y, x-y])) == "(x + y)^2 + (x - y)^2"
        @testset "SOSDecomposition equality" begin
            @polyvar x y
            @test !isapprox(SOSDecomposition([x+y, x-y]), SOSDecomposition([x+y]))
            @test !isapprox(SOSDecomposition([x+y, x-y]), SOSDecomposition([x+y, x+y]))
            @test isapprox(SOSDecomposition([x+y, x-y]), SOSDecomposition([x+y, x-y]))
            @test isapprox(SOSDecomposition([x+y, x-y]), SOSDecomposition([x-y, x+y+1e-8]), ztol=1e-7)
        end
    end
end
