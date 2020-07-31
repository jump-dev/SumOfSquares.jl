import MultivariatePolynomials
const MP = MultivariatePolynomials

import MultivariateBases
const MB = MultivariateBases

@testset "Monomial selection for certificate" begin
    @polyvar x y z
    @ncpolyvar a b
    @testset "Multipartite error not commutative" for parts in [([a],), ([a], [b])]
        err = ArgumentError("Multipartite Newton polytope not supported with noncommutative variables.")
        @test_throws err SumOfSquares.Certificate.monomials_half_newton_polytope([a*b, b^2], parts)
    end
    @testset "Multipartite error not disjoint: $parts" for parts in [([x], [x]),
                                                                     ([x], [x, y]),
                                                                     ([x], [y], [x, y])]
        err = ArgumentError("Parts are not disjoint in multipartite Newton polytope estimation: $parts.")
        @test_throws err SumOfSquares.Certificate.monomials_half_newton_polytope([x*y, y^2], parts)
    end
    @testset "Unipartite" begin
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([x*y, y^2], tuple()) == [y]
        @test isempty(SumOfSquares.Certificate.monomials_half_newton_polytope([x, y], tuple()))
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([x^2, y^2], tuple()) == [x, y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([x^2, y^2], ([x, y],)) == [x, y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([x^2, y^2], ([y, x],)) == [x, y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([x^2, x^3*y^2, x^4*y^4], tuple()) == [x^2*y^2, x]
    end
    @testset "Non-commutative" begin
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([a^4, a^3*b, a*b*a^2, a*b*a*b], tuple()) == [a^2, a*b]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([a^2, a^10*b^20*a^11, a^11*b^20*a^10, a^10*b^20*a^20*b^20*a^10], tuple()) == [a^10*b^20*a^10, a]
    end
    @testset "Multipartite" begin
        # In the part [y, z], the degree is between 0 and 2
        X = [x^4, x^2*y^2, x^2*z^2, x^2*y*z, y*z]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(X, tuple(), apply_post_filter = false) == [x^2, x*y, x*z, y*z, x, y, z]
        function full_test(X, Y, part1, part2)
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(X, (part1,), apply_post_filter = false) == Y
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(X, (part2,), apply_post_filter = false) == Y
            a = SumOfSquares.Certificate.monomials_half_newton_polytope(X, (part2,), apply_post_filter = false)
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(X, (part1, part2), apply_post_filter = false) == Y
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(X, (part2, part1), apply_post_filter = false) == Y
        end
        full_test(X, monovec([x^2, x*y, x*z,      x, y, z]), [x], [y, z])
        full_test(X, monovec([x^2, x*y, x*z, y*z, x,    z]), [y], [x, z])
        full_test(X, monovec([x^2, x*y, x*z, y*z, x, y   ]), [z], [x, y])
        # FIXME: With recursive merging, it should give [x^2, x*y, x*z, x]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope([x^4, x^2*y^2, x^2*z^2, x^2*y*z, y*z], ([x], [y], [z]), apply_post_filter = false) == [x^2, x*y, x*z, y*z, x, y, z]
    end
end

@testset "Random SOS should be SOS" begin
    @polyvar x y
    x = [1, x, y, x^2, y^2, x*y]
    @test_throws ArgumentError randsos(x, monotype=:Unknown)
    for i in 1:10
        for monotype in [:Classic, :Gram]
            p = randsos(x, monotype=monotype)
            @test p isa GramMatrix
            @test isposdef(Matrix(p.Q))
        end
    end
end

function certificate_api(certificate::Certificate.AbstractIdealCertificate)
    @polyvar x
    poly = x + 1
    domain = @set x == 1
    @test Certificate.get(certificate, Certificate.ReducedPolynomial(), poly, domain) isa MP.AbstractPolynomial
    basis = Certificate.get(certificate, Certificate.GramBasis(), poly)
    @test basis isa MB.AbstractPolynomialBasis || basis isa Vector{<:MB.AbstractPolynomialBasis}
    @test basis isa Certificate.get(typeof(certificate), Certificate.GramBasisType())
    @test Certificate.get(certificate, Certificate.Cone()) isa SumOfSquares.SOSLikeCone
    @test SumOfSquares.matrix_cone_type(typeof(certificate)) <: MOI.AbstractVectorSet
    zbasis = Certificate.zero_basis(certificate)
    @test zbasis <: MB.AbstractPolynomialBasis
    @test zbasis == Certificate.zero_basis_type(typeof(certificate))
end

@testset "API" begin
    @polyvar x
    cone = SumOfSquares.SOSCone()
    BT = MB.MonomialBasis
    basis = BT([x^2, x])
    @testset "$(typeof(certificate))" for certificate in [
        Certificate.MaxDegree(cone, BT, 2),
        Certificate.FixedBasis(cone, basis),
        Certificate.Newton(cone, BT, tuple())
    ]
        certificate_api(certificate)
        certificate_api(Certificate.Remainder(certificate))
        if certificate isa Certificate.MaxDegree
            certificate_api(Certificate.SparseIdeal(VariableSparsity(), certificate))
        end
        @testset "$(typeof(sparsity))" for sparsity in [SignSymmetry(), MonomialSparsity(1)]
            certificate_api(Certificate.SparseIdeal(sparsity, certificate))
            certificate_api(Certificate.SparseIdeal(sparsity, Certificate.Remainder(certificate)))
        end
    end
end

include("ceg_test.jl")
include("csp_test.jl")
include("sparsity.jl")
