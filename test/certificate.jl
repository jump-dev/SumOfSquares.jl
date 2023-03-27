import MultivariatePolynomials as MP

import MultivariateBases as MB

@testset "_merge_sorted" begin
    @test SumOfSquares.Certificate._merge_sorted([4, 1], [3, 0]) == [4, 3, 1, 0]
    @test SumOfSquares.Certificate._merge_sorted((4, 1), (3, 0)) == (4, 3, 1, 0)
    @test SumOfSquares.Certificate._merge_sorted([4, 1], [3, 2]) == [4, 3, 2, 1]
    @test SumOfSquares.Certificate._merge_sorted((4, 1), (3, 2)) == (4, 3, 2, 1)
end

@testset "with_variables" begin
    @polyvar x y z
    p = x + z
    v = SumOfSquares.Certificate.with_variables(p, y)
    @test v.inner === p
    @test MP.variables(v.inner) == [x, z]
    @test MP.variables(v) == [x, y, z]

    @ncpolyvar a b
    q = a * b + b * a
    v = SumOfSquares.Certificate.with_variables(q, FullSpace())
    @test v.inner === q
    @test MP.variables(v.inner) == [a, b, a]
    @test MP.variables(v) == [a, b, a]
end

@testset "Monomial selection for certificate" begin
    @polyvar x y z
    @ncpolyvar a b
    @testset "Multipartite error not commutative" for parts in
                                                      [([a],), ([a], [b])]
        err = ArgumentError(
            "Multipartite Newton polytope not supported with noncommutative variables.",
        )
        @test_throws err SumOfSquares.Certificate.monomials_half_newton_polytope(
            [a * b, b^2],
            Certificate.NewtonDegreeBounds(parts),
        )
    end
    @testset "Multipartite error not disjoint: $parts" for parts in [
        ([x], [x]),
        ([x], [x, y]),
        ([x], [y], [x, y]),
    ]
        err = ArgumentError(
            "Parts are not disjoint in multipartite Newton polytope estimation: $parts.",
        )
        @test_throws err SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x * y, y^2],
            Certificate.NewtonDegreeBounds(parts),
        )
    end
    uni = Certificate.NewtonDegreeBounds(tuple())
    @testset "Unipartite" begin
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x * y, y^2],
            uni,
        ) == [y]
        @test isempty(
            SumOfSquares.Certificate.monomials_half_newton_polytope(
                [x, y],
                uni,
            ),
        )
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x^2, y^2],
            uni,
        ) == [x, y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x^2, y^2],
            Certificate.NewtonDegreeBounds(([x, y],)),
        ) == [x, y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x^2, y^2],
            Certificate.NewtonDegreeBounds(([y, x],)),
        ) == [x, y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x^2, x^3 * y^2, x^4 * y^4],
            uni,
        ) == [x^2 * y^2, x, x * y, x^2, x * y^2, x^2 * y]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x^2, x^3 * y^2, x^4 * y^4],
            Certificate.NewtonFilter(uni),
        ) == [x^2 * y^2, x]
    end
    @testset "Non-commutative" begin
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [a^4, a^3 * b, a * b * a^2, a * b * a * b],
            uni,
        ) == [a^2, a * b]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [
                a^2,
                a^10 * b^20 * a^11,
                a^11 * b^20 * a^10,
                a^10 * b^20 * a^20 * b^20 * a^10,
            ],
            Certificate.NewtonFilter(uni),
        ) == [a^10 * b^20 * a^10, a]
    end
    @testset "Multipartite" begin
        # In the part [y, z], the degree is between 0 and 2
        X = [x^4, x^2 * y^2, x^2 * z^2, x^2 * y * z, y * z]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(X, uni) ==
              [x^2, x * y, x * z, y * z, x, y, z]
        function full_test(X, Y, part1, part2)
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part1,)),
            ) == Y
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part2,)),
            ) == Y
            a = SumOfSquares.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part2,)),
            )
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part1, part2)),
            ) == Y
            @test SumOfSquares.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part2, part1)),
            ) == Y
        end
        full_test(X, monovec([x^2, x * y, x * z, x, y, z]), [x], [y, z])
        full_test(X, monovec([x^2, x * y, x * z, y * z, x, z]), [y], [x, z])
        full_test(X, monovec([x^2, x * y, x * z, y * z, x, y]), [z], [x, y])
        # FIXME: With recursive merging, it should give [x^2, x*y, x*z, x]
        @test SumOfSquares.Certificate.monomials_half_newton_polytope(
            [x^4, x^2 * y^2, x^2 * z^2, x^2 * y * z, y * z],
            Certificate.NewtonDegreeBounds(([x], [y], [z])),
        ) == [x^2, x * y, x * z, y * z, x, y, z]
    end
end

@testset "Random SOS should be SOS" begin
    @polyvar x y
    x = [1, x, y, x^2, y^2, x * y]
    @test_throws ArgumentError randsos(x, monotype = :Unknown)
    for i in 1:10
        for monotype in [:Classic, :Gram]
            p = randsos(x, monotype = monotype)
            @test p isa GramMatrix
            @test isposdef(Matrix(p.Q))
        end
    end
end

function _certificate_api(certificate::Certificate.AbstractCertificate)
    @test Certificate.cone(certificate) isa SumOfSquares.SOSLikeCone
    @test SumOfSquares.matrix_cone_type(typeof(certificate)) <:
          MOI.AbstractVectorSet
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
        @test issorted(monos, rev = true)
    end
end
function _basis_check(basis, basis_type)
    @test basis isa MB.AbstractPolynomialBasis ||
          basis isa Vector{<:MB.AbstractPolynomialBasis}
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
    @test Certificate.reduced_polynomial(certificate, poly, domain) isa
          MP.AbstractPolynomial
    _basis_check(
        Certificate.gram_basis(certificate, poly),
        Certificate.gram_basis_type(typeof(certificate)),
    )
    zbasis = Certificate.zero_basis(certificate)
    @test zbasis <: MB.AbstractPolynomialBasis
    @test zbasis == Certificate.zero_basis_type(typeof(certificate))
end

function certificate_api(certificate::Certificate.AbstractPreorderCertificate)
    _certificate_api(certificate)
    @polyvar x
    poly = x + 1
    domain = @set x >= 1
    processed = Certificate.preprocessed_domain(certificate, domain, poly)
    for idx in Certificate.preorder_indices(certificate, processed)
        _basis_check(
            Certificate.multiplier_basis(certificate, idx, processed),
            Certificate.multiplier_basis_type(typeof(certificate)),
        )
        @test Certificate.generator(certificate, idx, processed) isa
              MP.AbstractPolynomial
    end
    icert = Certificate.ideal_certificate(certificate)
    @test icert isa Certificate.AbstractIdealCertificate
    @test typeof(icert) == Certificate.ideal_certificate(typeof(certificate))
end

@testset "API" begin
    @polyvar x
    cone = SumOfSquares.SOSCone()
    BT = MB.MonomialBasis
    maxdegree = 2
    function _test(certificate::Certificate.AbstractIdealCertificate)
        certificate_api(certificate)
        mult_cert = certificate
        if mult_cert isa Certificate.Sparsity.Ideal
            mult_cert = mult_cert.certificate
        end
        if mult_cert isa Certificate.Remainder # FIXME not supported yet as mult cert
            mult_cert = mult_cert.gram_certificate
        end
        if mult_cert isa Certificate.FixedBasis # FIXME not supported yet
            mult_cert = Certificate.MaxDegree(cone, BT, maxdegree)
        end
        preorder = Certificate.Putinar(mult_cert, certificate, maxdegree)
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
        Certificate.Newton(cone, BT, tuple()),
    ]
        _test(certificate)
        _test(Certificate.Remainder(certificate))
        if certificate isa Certificate.MaxDegree
            _test(Certificate.Sparsity.Ideal(Sparsity.Variable(), certificate))
        end
        @testset "$(typeof(sparsity))" for sparsity in [
            SignSymmetry(),
            Sparsity.Monomial(ChordalCompletion(), 1),
        ]
            _test(Certificate.Sparsity.Ideal(sparsity, certificate))
            _test(
                Certificate.Sparsity.Ideal(
                    sparsity,
                    Certificate.Remainder(certificate),
                ),
            )
        end
    end
end

function test_putinar_ijk(i, j, k, default::Bool, post_filter::Bool = default)
    v = @polyvar x y
    poly = x^(2i) + y^(2j + 1)
    domain = @set y^(2k + 1) >= 0
    if default
        certificate =
            JuMP.moi_set(SOSCone(), monomials(poly); domain).certificate
    else
        newton = Certificate.NewtonDegreeBounds(tuple())
        if post_filter
            newton = Certificate.NewtonFilter(newton)
        end
        cert = Certificate.Newton(SOSCone(), MB.MonomialBasis, newton)
        certificate = Certificate.Putinar(cert, cert, max(2i, 2j + 1, 2k + 1))
    end
    processed = Certificate.preprocessed_domain(certificate, domain, poly)
    for idx in Certificate.preorder_indices(certificate, processed)
        monos =
            Certificate.multiplier_basis(certificate, idx, processed).monomials
        if k > j
            @test isempty(monos)
        else
            w = post_filter ? v[2:2] : v
            @test monos == MP.monomials(w, max(0, min(i, j) - k):(j-k))
        end
    end
    icert = Certificate.ideal_certificate(certificate)
    @test icert isa Certificate.Newton
end

@testset "Putinar $i $j $k" for (i, j, k) in [(1, 1, 2), (1, 3, 2), (3, 2, 1)] #, (4, 2, 1)]
    @testset "post_filter=$post" for post in [true, false]
        test_putinar_ijk(i, j, k, post)
    end
end

include("ceg_test.jl")
include("csp_test.jl")
include("sparsity.jl")
include("symmetry.jl")
