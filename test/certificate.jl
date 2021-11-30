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

function _certificate_api(certificate::Certificate.AbstractCertificate)
    @test Certificate.get(certificate, Certificate.Cone()) isa SumOfSquares.SOSLikeCone
    @test SumOfSquares.matrix_cone_type(typeof(certificate)) <: MOI.AbstractVectorSet
end
function _basis_check_each(basis::MB.AbstractPolynomialBasis, basis_type)
    @test basis isa basis_type
    if basis isa MB.AbstractMonomialBasis
        # This fails if `basis` is `Vector{Monomial{true}}` instead of `MonomialVector{true}`
        # for DynamicPolynomials. This is important as
        # `polynomial(::AbstractMatrix, ::MonomialVector, ::Type)` is implemented but
        # `polynomial(Q::AbstractMatrix, X::AbstractVector, ::Type)` falls back to
        # `dot(X, Q * X)` which does not work (`promote_operation` calls `zero(eltype(X))`
        # which gives `Polynomial{true, Int}` which then tries to multiply a
        # `ScalarAffineFunction{Float64}` with an `Int`).
        monos = basis.monomials
        @test typeof(monos) == typeof(monovec(monos))
        @test issorted(monos, rev=true)
    end
end
function _basis_check(basis, basis_type)
    @test basis isa MB.AbstractPolynomialBasis || basis isa Vector{<:MB.AbstractPolynomialBasis}
    if basis isa Vector
        # FIXME `basis_type` is `Vector{MB.MonomialBasis}` instead of `Vector{MB.MonomialBasis{...}}`
        # Once this is fixed, we should check
        # @test basis isa basis_type
        for b in basis
            _basis_check_each(b, eltype(basis_type))
        end
    else
        _basis_check_each(basis, basis_type)
    end
end

function certificate_api(certificate::Certificate.AbstractIdealCertificate)
    _certificate_api(certificate)
    @polyvar x
    poly = x + 1
    domain = @set x == 1
    @test Certificate.get(certificate, Certificate.ReducedPolynomial(), poly, domain) isa MP.AbstractPolynomial
    _basis_check(Certificate.get(certificate, Certificate.GramBasis(), poly),
                 Certificate.get(typeof(certificate), Certificate.GramBasisType()))
    zbasis = Certificate.zero_basis(certificate)
    @test zbasis <: MB.AbstractPolynomialBasis
    @test zbasis == Certificate.zero_basis_type(typeof(certificate))
end

function certificate_api(certificate::Certificate.AbstractPreorderCertificate)
    _certificate_api(certificate)
    @polyvar x
    poly = x + 1
    domain = @set x >= 1
    processed = Certificate.get(certificate, Certificate.PreprocessedDomain(), domain, poly)
    for idx in Certificate.get(certificate, Certificate.PreorderIndices(), processed)
        _basis_check(Certificate.get(certificate, Certificate.MultiplierBasis(), idx, processed),
                     Certificate.get(typeof(certificate), Certificate.MultiplierBasisType()))
        @test Certificate.get(certificate, Certificate.Generator(), idx, processed) isa MP.AbstractPolynomial
    end
    icert = Certificate.get(certificate, Certificate.IdealCertificate())
    @test icert isa Certificate.AbstractIdealCertificate
    @test typeof(icert) == Certificate.get(typeof(certificate), Certificate.IdealCertificate())
end


@testset "API" begin
    @polyvar x
    cone = SumOfSquares.SOSCone()
    BT = MB.MonomialBasis
    maxdegree = 2
    function _test(certificate::Certificate.AbstractIdealCertificate)
        certificate_api(certificate)
        preorder = Certificate.Putinar(certificate, cone, BT, maxdegree)
        certificate_api(preorder)
        sparsities = Sparsity.Pattern[Sparsity.Variable()]
        if certificate isa Certificate.MaxDegree
            push!(sparsities, Sparsity.Monomial(ChordalCompletion(), 1))
        end
        @testset "$(typeof(sparsity))" for sparsity in sparsities
            certificate_api(Certificate.Sparsity.Preorder(sparsity, preorder))
        end
    end
    basis = BT(monovec([x^2, x]))
    @testset "$(typeof(certificate))" for certificate in [
        Certificate.MaxDegree(cone, BT, maxdegree),
        Certificate.FixedBasis(cone, basis),
        Certificate.Newton(cone, BT, tuple())
    ]
        _test(certificate)
        _test(Certificate.Remainder(certificate))
        if certificate isa Certificate.MaxDegree
            _test(Certificate.Sparsity.Ideal(Sparsity.Variable(), certificate))
        end
        @testset "$(typeof(sparsity))" for sparsity in [SignSymmetry(), Sparsity.Monomial(ChordalCompletion(), 1)]
            _test(Certificate.Sparsity.Ideal(sparsity, certificate))
            _test(Certificate.Sparsity.Ideal(sparsity, Certificate.Remainder(certificate)))
        end
    end
end

include("ceg_test.jl")
include("csp_test.jl")
include("sparsity.jl")
