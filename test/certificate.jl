import StarAlgebras as SA
import MutableArithmetics as MA
import MultivariatePolynomials as MP
import MultivariateBases as MB

const SOS = SumOfSquares

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
        @test_throws err SOS.Certificate.monomials_half_newton_polytope(
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
        @test_throws err SOS.Certificate.monomials_half_newton_polytope(
            [x * y, y^2],
            Certificate.NewtonDegreeBounds(parts),
        )
    end
    uni = Certificate.NewtonDegreeBounds(tuple())
    @testset "Unipartite" begin
        @test SOS.Certificate.monomials_half_newton_polytope(
            [x * y, y^2],
            uni,
        ) == [y]
        @test isempty(
            SOS.Certificate.monomials_half_newton_polytope([x, y], uni),
        )
        @test SOS.Certificate.monomials_half_newton_polytope([x^2, y^2], uni) ==
              [x, y]
        @test SOS.Certificate.monomials_half_newton_polytope(
            [x^2, y^2],
            Certificate.NewtonDegreeBounds(([x, y],)),
        ) == [x, y]
        @test SOS.Certificate.monomials_half_newton_polytope(
            [x^2, y^2],
            Certificate.NewtonDegreeBounds(([y, x],)),
        ) == [x, y]
        @test SOS.Certificate.monomials_half_newton_polytope(
            [x^2, x^3 * y^2, x^4 * y^4],
            uni,
        ) == [x^2 * y^2, x, x * y, x^2, x * y^2, x^2 * y]
        @test SOS.Certificate.monomials_half_newton_polytope(
            [x^2, x^3 * y^2, x^4 * y^4],
            Certificate.NewtonFilter(uni),
        ) == [x^2 * y^2, x]
    end
    @testset "Non-commutative" begin
        @test SOS.Certificate.monomials_half_newton_polytope(
            [a^4, a^3 * b, a * b * a^2, a * b * a * b],
            uni,
        ) == [a^2, a * b]
        @test SOS.Certificate.monomials_half_newton_polytope(
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
        @test SOS.Certificate.monomials_half_newton_polytope(X, uni) ==
              [x^2, x * y, x * z, y * z, x, y, z]
        function full_test(X, Y, part1, part2)
            @test SOS.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part1,)),
            ) == Y
            @test SOS.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part2,)),
            ) == Y
            a = SOS.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part2,)),
            )
            @test SOS.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part1, part2)),
            ) == Y
            @test SOS.Certificate.monomials_half_newton_polytope(
                X,
                Certificate.NewtonDegreeBounds((part2, part1)),
            ) == Y
        end
        full_test(X, monomial_vector([x^2, x * y, x * z, x, y, z]), [x], [y, z])
        full_test(
            X,
            monomial_vector([x^2, x * y, x * z, y * z, x, z]),
            [y],
            [x, z],
        )
        full_test(
            X,
            monomial_vector([x^2, x * y, x * z, y * z, x, y]),
            [z],
            [x, y],
        )
        # FIXME: With recursive merging, it should give [x^2, x*y, x*z, x]
        @test SOS.Certificate.monomials_half_newton_polytope(
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
function _basis_check_each(basis::SA.ExplicitBasis, basis_type)
    if basis isa MB.SubBasis
        # This fails if `basis` is `Vector{<:Monomial}` instead of `MonomialVector`
        # for DynamicPolynomials. This is important as
        # `polynomial(::AbstractMatrix, ::MonomialVector, ::Type)` is implemented but
        # `polynomial(Q::AbstractMatrix, X::AbstractVector, ::Type)` falls back to
        # `dot(X, Q * X)` which does not work (`promote_operation` calls `zero(eltype(X))`
        # which gives `Polynomial{true, Int}` which then tries to multiply a
        # `ScalarAffineFunction{Float64}` with an `Int`).
        monos = basis.monomials
        @test typeof(monos) == typeof(monomial_vector(monos))
        @test issorted(monos)
    end
end
function _basis_check(basis, basis_type)
    @test basis isa SA.ExplicitBasis || basis isa Vector{<:SA.ExplicitBasis}
    @test basis isa basis_type
    if basis isa Vector
        for b in basis
            _basis_check_each(b, basis_type)
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
    a = MB.algebra_element(
        MB.sparse_coefficients(poly),
        MB.FullBasis{MB.Monomial,MP.monomial_type(poly)}(),
    )
    @test Certificate.reduced_polynomial(certificate, a, domain) isa
          SA.AlgebraElement
    gram_basis = Certificate.gram_basis(
        certificate,
        Certificate.WithVariables(a, MP.variables(poly)),
    )
    _basis_check(
        gram_basis,
        MA.promote_operation(Certificate.gram_basis, typeof(certificate)),
    )
    flat_bases, flat_weights, _ = SumOfSquares.Bridges.Constraint._flatten([gram_basis], [a])
    zbasis = Certificate.zero_basis(certificate, MB.explicit_basis(a), domain, flat_bases, flat_weights)
    @test zbasis isa SA.AbstractBasis
    @test zbasis isa MA.promote_operation(Certificate.zero_basis, typeof(certificate), typeof(MB.explicit_basis(a)), typeof(domain), typeof(gram_basis), Vector{typeof(a)})
end

function certificate_api(certificate::Certificate.AbstractPreorderCertificate)
    _certificate_api(certificate)
    @polyvar x
    poly = x + 1
    domain = @set x >= 1
    a = MB.algebra_element(
        MB.sparse_coefficients(poly),
        MB.FullBasis{MB.Monomial,MP.monomial_type(poly)}(),
    )
    processed = Certificate.preprocessed_domain(
        certificate,
        domain,
        Certificate.WithVariables(a, MP.variables(poly)),
    )
    for idx in Certificate.preorder_indices(certificate, processed)
        _basis_check(
            Certificate.multiplier_basis(certificate, idx, processed),
            Certificate.multiplier_basis_type(
                typeof(certificate),
                MP.monomial_type(x),
            ),
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
    B = MB.Monomial
    full_basis = MB.FullBasis{B,MP.monomial_type(x)}()
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
            mult_cert = Certificate.MaxDegree(cone, full_basis, full_basis, maxdegree)
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
    basis = MB.SubBasis{B}(monomial_vector([x^2, x]))
    @testset "$(nameof(typeof(certificate)))" for certificate in [
        Certificate.MaxDegree(cone, full_basis, full_basis, maxdegree),
        Certificate.FixedBasis(cone, basis, full_basis),
        Certificate.Newton(cone, full_basis, full_basis, tuple()),
    ]
        _test(certificate)
        _test(Certificate.Remainder(certificate))
        if certificate isa Certificate.MaxDegree
            _test(Certificate.Sparsity.Ideal(Sparsity.Variable(), certificate))
        end
        @testset "$(nameof(typeof(sparsity)))" for sparsity in [
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
            JuMP.moi_set(
                SOSCone(),
                MB.SubBasis{MB.Monomial}(monomials(poly)),
                MB.FullBasis{MB.Monomial,MP.monomial_type(poly)}(),
                MB.FullBasis{MB.Monomial,MP.monomial_type(poly)}();
                domain,
            ).certificate
    else
        newton = Certificate.NewtonDegreeBounds(tuple())
        if post_filter
            newton = Certificate.NewtonFilter(newton)
        end
        cert = Certificate.Newton(
            SOSCone(),
            MB.FullBasis{MB.Monomial,MP.monomial_type(x * y)}(),
            MB.FullBasis{MB.Monomial,MP.monomial_type(x * y)}(),
            newton,
        )
        certificate = Certificate.Putinar(cert, cert, max(2i, 2j + 1, 2k + 1))
    end
    alg_el = MB.algebra_element(
        MB.sparse_coefficients(poly),
        MB.FullBasis{MB.Monomial,MP.monomial_type(poly)}(),
    )
    processed = Certificate.preprocessed_domain(certificate, domain, alg_el)
    for idx in Certificate.preorder_indices(certificate, processed)
        monos =
            Certificate.multiplier_basis(certificate, idx, processed).monomials
        if k > j
            @test isempty(monos)
        else
            w = post_filter ? v[2:2] : v
            @test monos == MP.monomials(
                w,
                max(0, (post_filter ? j : min(i, j)) - k):(j-k),
            )
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
