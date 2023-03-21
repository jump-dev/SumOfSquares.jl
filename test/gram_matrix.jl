using LinearAlgebra, Test, SumOfSquares

@testset "GramMatrix tests" begin
    @testset "GramMatrix" begin
        @polyvar x y
        P = GramMatrix{Int}((i, j) -> i + j, [x^2, x * y, y^2])
        @test nvariables(P) == 2
        @test variables(P) == [x, y]
        Pp1 = 2x^4 + 6x^3 * y + 12x^2 * y^2 + 10x * y^3 + 6y^4 + 1
        @test P + 1 == Pp1
        @test P + 1.0 == Pp1
        zP = zero(typeof(P))
        @test isempty(zP.Q)
        @test zP == 0
        p = polynomial(P)
        @test coefficients(p) == [2, 6, 12, 10, 6]
        @test monomialtype(p) == typeof(x * y)
        @test monomialtype(typeof(p)) == typeof(x * y)
        @test monomials(p) == [x^4, x^3 * y, x^2 * y^2, x * y^3, y^4]
        for i in 1:3
            for j in 1:3
                @test P[i, j] == i + j
            end
        end
        for P in (
            GramMatrix{Int}((i, j) -> i * j, [y, x]),
            GramMatrix{Int}((i, j) -> (3 - i) * (3 - j), monovec([y, x])),
            GramMatrix([1 2; 2 4], [y, x]),
            GramMatrix([4 2; 2 1], monovec([y, x])),
        )
            @test P.Q.Q == [4, 2, 1]
            @test P.basis.monomials[1] == x
            @test P.basis.monomials[2] == y
        end
        P = GramMatrix{Int}(
            (i, j) -> ((i, j) == (1, 1) ? 2 : 0),
            [x * y, x^2, y^2],
        )
        Q = GramMatrix([0 1; 1 0], [x^2, y^2])
        @test P == Q
        p = GramMatrix([2 3; 3 2], [x, y])
        @test polynomial(p) isa AbstractPolynomial{Int}
        @test polynomial(p, Float64) isa AbstractPolynomial{Float64}

        @test -2 * x +
              dot(-x - x^2, 0) +
              GramMatrix{Int}((i, j) -> 1, [1, x]) == -(-x^2 - 1)
        P = GramMatrix{Int}((i, j) -> i + j, [x^2, x * y, y^2])
        p = polynomial(P)
        @test !iszero(P)
        @test iszero(P - P)
        @test iszero(P - p)
        @test iszero(p - P)
        @test P + P == P + p
        @test x * P == x * p
        @test 2 * P == 2 * p
        @test P * (2x) == (2x) * p

        @test differentiate(GramMatrix{Int}((i, j) -> 1, [x]), y) == 0
        @test differentiate(GramMatrix{Int}((i, j) -> 1, [x, y]), y, 1) ==
              2x + 2y
        #@inferred differentiate(GramMatrix{Int}((i,j)->1, [x, y]), y, 0) # FIXME failing at the moment

        @test GramMatrix([2 3; 3 2], [x, y]) == 2x^2 + 2y^2 + 6x * y

        mp = GramMatrix{Int}((i, j) -> i + j, [x])
        @test norm(mp) == norm(mp, 2) == 2.0

        @polyvar v[1:3]
        P = [1 2 3; 2 4 5; 3 5 6]
        p = GramMatrix(P, v)
        #@inferred p(x => ones(3))
        @test p(v => ones(3)) == 31
        #@inferred subs(p, x => ones(3))
        @test subs(p, v => ones(3)) == 31
        @testset "Gram Operate" begin
            p = GramMatrix([2 3; 3 2], [x, y])
            q = GramMatrix(5 * ones(1, 1), [x])
            r = @inferred gram_operate(/, q, 5)
            @test r.Q == ones(1, 1)
            @test r.basis.monomials == [x]
            r = @inferred gram_operate(+, p, q)
            @test r.Q == [7 3; 3 2]
            @test r.basis.monomials == [x, y]
            q = GramMatrix(5 * ones(1, 1), [y])
            r = @inferred gram_operate(+, p, q)
            @test r.Q == [2 3; 3 7]
            @test r.basis.monomials == [x, y]
            q = GramMatrix([5.0 7; 7 9], [x * y, 1])
            r = @inferred gram_operate(+, p, q)
            @test r.Q == [
                5 0 0 7
                0 2 3 0
                0 3 2 0
                7 0 0 9
            ]
            @test r.basis.monomials == [x * y, x, y, 1]
        end
        @testset "With VariableIndex" begin
            a = MOI.VariableIndex(1)
            g = GramMatrix(SymMatrix([a, a, a], 2), [x, y])
            U = MOI.ScalarAffineFunction{Float64}
            @test coefficienttype(g) == U
            @test g isa AbstractPolynomialLike{U}
            @test polynomialtype(g) <: AbstractPolynomial{U}
            @test polynomialtype(typeof(g)) <: AbstractPolynomial{U}
            @test polynomial(g) isa AbstractPolynomial{U}
        end
    end
    @testset "SOSDecomposition" begin
        @polyvar x y
        #       @test isempty(SOSDecomposition(typeof(x)[]))
        #       ps = [1, x + y, x^2, x*y, 1 + x + x^2]
        #       P = GramMatrix(SOSDecomposition(ps))
        #       P.Q == [2 0 1 0 1; 0 1 0 0 0; 1 0 2 1 1; 0 0 1 1 0; 1 0 1 0 2]
        #       P.basis.monomials == [x^2, x*y, x, y, 1]
        #       @test P == P
        #       @test isapprox(GramMatrix(SOSDecomposition(P)), P)
        P = GramMatrix{Int}((i, j) -> i + j, [x^2, x * y, y^2])
        @test polynomialtype(SOSDecomposition(P)) <: AbstractPolynomialLike
        @test sprint(show, SOSDecomposition([x + y, x - y])) ==
              "(x + y)^2 + (x - y)^2"
        @test polynomial(SOSDecomposition([x + y, x - y])) ==
              (x + y)^2 + (x - y)^2
        @test polynomial(SOSDecomposition([x + y, x - y]), Float64) ==
              (x + y)^2 + (x - y)^2
        @testset "SOSDecomposition equality" begin
            @polyvar x y
            @test !isapprox(
                SOSDecomposition([x + y, x - y]),
                SOSDecomposition([x + y]),
            )
            @test !isapprox(
                SOSDecomposition([x + y, x - y]),
                SOSDecomposition([x + y, x + y]),
            )
            @test isapprox(
                SOSDecomposition([x + y, x - y]),
                SOSDecomposition([x + y, x - y]),
            )
            @test isapprox(
                SOSDecomposition([x + y, x - y]),
                SOSDecomposition([x - y, x + y + 1e-8]),
                ztol = 1e-7,
            )
        end
        @testset "With VariableIndex" begin
            a = MOI.VariableIndex(1)
            p = polynomial([a], [x])
            q = polynomial([a], [y])
            s = SOSDecomposition([p, q])
            U = MOI.ScalarQuadraticFunction{Float64}
            @test coefficienttype(s) == U
            @test s isa AbstractPolynomialLike{U}
            @test polynomialtype(s) <: AbstractPolynomial{U}
            @test polynomialtype(typeof(s)) <: AbstractPolynomial{U}
        end
    end

    @testset "SOSDecompositionWithDomain" begin
        @polyvar x y
        K = @set 1 - x^2 >= 0 && 1 - y^2 >= 0
        ps = SOSDecomposition([x + y, x - y])
        ps1 = SOSDecomposition([x])
        ps2 = SOSDecomposition([y])
        @test [ps, ps1] isa Vector{
            SOSDecomposition{Int,T,Int},
        } where {T<:AbstractPolynomialLike}
        @test sprint(show, SOSDecompositionWithDomain(ps, [ps1, ps2], K)) ==
              "(x + y)^2 + (x - y)^2 + (x)^2 * (-x^2 + 1) + (y)^2 * (-y^2 + 1)"

        @testset "SOSDecompositionWithDomain equality" begin
            @polyvar x y
            K = @set 1 - x^2 >= 0 && 1 - y^2 >= 0
            B = @set 1 - x >= 0 && 1 - y >= 0
            ps = SOSDecomposition([x + y, x - y])
            ps1 = SOSDecomposition([x + y, x^2 - y])
            ps2 = SOSDecomposition([x + y, y^2 - x])
            sosdec = SOSDecompositionWithDomain(ps, [ps1, ps2], K)
            @test typeof(polynomial(sosdec)) <: AbstractPolynomialLike
            @test isapprox(sosdec, sosdec)
            @test !isapprox(
                sosdec,
                SOSDecompositionWithDomain(ps, [ps1, ps2], B),
            )
            @test !isapprox(
                SOSDecompositionWithDomain(ps, [ps1, ps1], K),
                sosdec,
            )
        end
    end

    @testset "build_gram_matrix" begin
        v = MOI.VariableIndex.(1:3)
        w = MOI.VariableIndex.(1:4)
        @polyvar x y
        basis = MonomialBasis(monomials([x, y], 1))
        @testset "$T" for T in [Float64, Int, BigFloat]
            #@test_throws DimensionMismatch SumOfSquares.build_gram_matrix(w, basis, T, MOI.PositiveSemidefiniteConeTriangle)
            g = SumOfSquares.build_gram_matrix(
                v,
                basis,
                MOI.PositiveSemidefiniteConeTriangle,
                T,
            )
            @test g isa GramMatrix{
                MOI.VariableIndex,
                typeof(basis),
                MOI.ScalarAffineFunction{T},
            }
            @test g.Q[1, 2] == v[2]
            p = polynomial(g)
            @test p isa AbstractPolynomial{MOI.ScalarAffineFunction{T}}
            @test typeof(p) == polynomialtype(g)
            #@test_throws DimensionMismatch SumOfSquares.build_gram_matrix(v, basis, T, SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle)
            h = SumOfSquares.build_gram_matrix(
                w,
                basis,
                SumOfSquares.COI.HermitianPositiveSemidefiniteConeTriangle,
                T,
            )
            @test h isa GramMatrix{
                MOI.ScalarAffineFunction{Complex{T}},
                typeof(basis),
                MOI.ScalarAffineFunction{Complex{T}},
                SumOfSquares.MultivariateMoments.VectorizedHermitianMatrix{
                    MOI.VariableIndex,
                    T,
                    MOI.ScalarAffineFunction{Complex{T}},
                },
            }
            @test h.Q[1, 2] ≈ w[2] + im * w[4]
            q = polynomial(h)
            @test q isa AbstractPolynomial{MOI.ScalarAffineFunction{Complex{T}}}
            @test typeof(q) == polynomialtype(h)
        end
    end
end
