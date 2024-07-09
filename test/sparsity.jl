using Test
import Combinatorics
using SumOfSquares
import MultivariateBases as MB

function _algebra_element(monos)
    return MB.algebra_element(
        SA.SparseCoefficients(monos, ones(length(monos))),
        MB.FullBasis{MB.Monomial,eltype(monos)}(),
    )
end

function _xor_complement_test(
    exps,
    expected,
    n = maximum(ndigits(exp, base = 2) for exp in exps; init = 1),
)
    for e in Combinatorics.permutations(exps)
        @test Certificate.Sparsity.xor_complement(e, n) == expected
    end
end

function xor_complement_test()
    _xor_complement_test([1], Int[])
    _xor_complement_test(Int[], [1])
    _xor_complement_test([1], [2], 2)
    _xor_complement_test([2], [1])
    _xor_complement_test([1, 2], Int[])
    _xor_complement_test([1, 3], Int[])
    _xor_complement_test(Int[], [1, 2], 2)
    _xor_complement_test([7], [3, 5])
    _xor_complement_test([5, 6, 3], [7])
    _xor_complement_test([3], [3, 4], 3)
    _xor_complement_test([32, 3, 24, 14, 21, 56], [7, 27])
    return
end

function set_monos(bases::Vector{<:MB.SubBasis})
    return Set([basis.monomials for basis in bases])
end

"""
    wml19()

Examples of [MWL19].

[WML19] Wang, Jie, Victor Magron, and Jean-Bernard Lasserre. "TSSOS: A Moment-SOS hierarchy that exploits term sparsity." arXiv preprint arXiv:1912.08899 (2019).
"""
function wml19()
    @polyvar x[1:3]
    certificate = Certificate.Newton(
        SOSCone(),
        MB.FullBasis{MB.Monomial,typeof(prod(x))}(),
        MB.FullBasis{MB.Monomial,typeof(prod(x))}(),
        tuple(),
    )
    @testset "Example 4.2" begin
        f = 1 + x[1]^4 + x[2]^4 + x[3]^4 + prod(x) + x[2]
        alg_el = MB.algebra_element(
            MB.sparse_coefficients(f),
            MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
        )
        with_var = SumOfSquares.Certificate.WithVariables(alg_el, x)
        expected_1_false = Set(
            monomial_vector.([
                [x[3]^2],
                [x[1] * x[3], x[2]],
                [x[2], 1],
                [x[2]^2],
                [x[2] * x[3], x[1]],
                [x[1] * x[2], x[3]],
                [x[1]^2],
            ]),
        )
        expected_2 = Set(
            monomial_vector.([
                [x[1]^2, x[2]^2, x[3]^2, 1],
                [x[2], 1],
                [x[2] * x[3], x[1]],
                [x[1] * x[3], x[2]],
                [x[1] * x[2], x[3]],
            ]),
        )
        cluster_expected_1_false = Set(
            monomial_vector.([
                [1, x[2], x[1] * x[3]],
                [x[1]^2],
                [x[2]^2],
                [x[3]^2],
                [x[1], x[2] * x[3]],
                [x[3], x[1] * x[2]],
            ]),
        )
        cluster_expected_1 = Set(
            monomial_vector.([
                [1, x[2], x[1]^2, x[2]^2, x[1] * x[3], x[3]^2],
                [x[1], x[2] * x[3]],
                [x[3], x[1] * x[2]],
            ]),
        )
        cluster_expected_2 = Set(
            monomial_vector.([
                [1, x[2], x[1]^2, x[2]^2, x[1] * x[3], x[3]^2],
                [x[1], x[2] * x[3], x[3], x[1] * x[2]],
            ]),
        )
        @testset "$(nameof(typeof(completion))) $k $use_all_monomials" for completion in
                                                                           [
                ClusterCompletion(),
                ChordalCompletion(),
            ],
            k in 0:2,
            use_all_monomials in [false, true]

            if completion isa ClusterCompletion
                expected =
                    k == 1 ?
                    (
                        use_all_monomials ? cluster_expected_1 :
                        cluster_expected_1_false
                    ) : cluster_expected_2
            else
                expected =
                    (k == 1 && !use_all_monomials) ? expected_1_false :
                    expected_2
            end
            @test set_monos(
                Certificate.Sparsity.sparsity(
                    with_var,
                    Sparsity.Monomial(completion, k, use_all_monomials),
                    certificate,
                ),
            ) == expected
        end
        expected = Set(
            monomial_vector.([
                [x[1]^2, x[1] * x[3], x[2]^2, x[3]^2, x[2], 1],
                [x[1] * x[2], x[2] * x[3], x[1], x[3]],
            ]),
        )
        @test set_monos(
            Certificate.Sparsity.sparsity(
                with_var,
                SignSymmetry(),
                certificate,
            ),
        ) == expected
    end
    @testset "Example 5.4 $(nameof(typeof(ideal_certificate)))" for ideal_certificate in
                                                                    [
        Certificate.MaxDegree(
            SOSCone(),
            MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
            MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
            4,
        ),
        Certificate.Newton(
            SOSCone(),
            MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
            MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
            tuple(),
        ),
    ]
        preorder_certificate = Certificate.Putinar(
            Certificate.MaxDegree(
                SOSCone(),
                MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
                MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
                4,
            ),
            ideal_certificate,
            4,
        )
        f = x[1]^4 + x[2]^4 + x[1] * x[2]
        alg_el = MB.algebra_element(
            MB.sparse_coefficients(f),
            MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
        )
        K = @set 1 - 2x[1]^2 - x[2]^2 >= 0
        @testset "$(nameof(typeof(completion))) $k $use_all_monomials" for completion in
                                                                           [
                ClusterCompletion(),
                ChordalCompletion(),
            ],
            k in 0:2,
            use_all_monomials in [false, true]

            basis, preorder_bases = Certificate.Sparsity.sparsity(
                alg_el,
                K,
                Sparsity.Monomial(completion, k, use_all_monomials),
                preorder_certificate,
            )
            if k == 1 &&
               (!use_all_monomials || completion isa ChordalCompletion)
                if use_all_monomials
                    @test set_monos(preorder_bases[1]) == Set(
                        monomial_vector.([
                            [x[1], x[2]],
                            [constant_monomial(x[1] * x[2])],
                        ]),
                    )
                    @test set_monos(basis) == Set(
                        monomial_vector.([
                            [x[1], x[2]],
                            [x[1] * x[2], 1],
                            [x[1]^2, x[2]^2, 1],
                        ]),
                    )
                else
                    @test set_monos(preorder_bases[1]) ==
                          Set(monomial_vector.([[x[1], x[2]]]))
                    @test set_monos(basis) == Set(
                        monomial_vector.([
                            [x[1]^2],
                            [x[1] * x[2], 1],
                            [x[2]^2],
                            [x[1], x[2]],
                        ]),
                    )
                end
            else
                @test set_monos(preorder_bases[1]) == Set(
                    monomial_vector.([
                        [constant_monomial(x[1] * x[2])],
                        [x[1], x[2]],
                    ]),
                )
                @test set_monos(basis) == Set(
                    monomial_vector.([
                        [x[1]^2, x[1] * x[2], x[2]^2, 1],
                        [x[1], x[2]],
                    ]),
                )
            end
        end
    end
    @testset "Example 6.7" begin
        f =
            1 + x[1]^2 * x[2]^4 + x[1]^4 * x[2]^2 + x[1]^4 * x[2]^4 -
            x[1] * x[2]^2 - 3x[1]^2 * x[2]^2
        alg_el = MB.algebra_element(
            MB.sparse_coefficients(f),
            MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
        )
        with_var = SumOfSquares.Certificate.WithVariables(alg_el, x)
        basis = MB.SubBasis{MB.Monomial}(MP.monomials(f))
        @testset "$(nameof(typeof(completion))) $k $use_all_monomials" for completion in
                                                                           [
                ClusterCompletion(),
                ChordalCompletion(),
            ],
            k in 0:2,
            use_all_monomials in [false, true]

            expected = if completion isa ClusterCompletion
                Set(
                    monomial_vector.([
                        [x[1]^2 * x[2]^2, x[1] * x[2]^2, 1],
                        [x[1] * x[2]],
                        [x[1]^2 * x[2]],
                    ]),
                )
            else
                Set(
                    monomial_vector.([
                        [x[1] * x[2]^2, 1],
                        [x[1]^2 * x[2]^2, 1],
                        [x[1] * x[2]],
                        [x[1]^2 * x[2]],
                    ]),
                )
            end
            @test set_monos(
                Certificate.Sparsity.sparsity(
                    with_var,
                    Sparsity.Monomial(completion, k, use_all_monomials),
                    certificate,
                ),
            ) == expected
        end
        @test set_monos(
            Certificate.Sparsity.sparsity(
                with_var,
                SignSymmetry(),
                certificate,
            ),
        ) == Set(
            monomial_vector.([
                [x[1]^2 * x[2]^2, x[1] * x[2]^2, 1, x[1]^2 * x[2], x[1] * x[2]],
            ]),
        )
    end
end

"""
    l09()

Examples of [MWL19].

[L09] Lofberg, Johan. "Pre-and post-processing sum-of-squares programs in practice." IEEE transactions on automatic control 54.5 (2009): 1007-1011.
"""
function l09()
    @polyvar x[1:3]
    certificate = Certificate.Newton(
        SOSCone(),
        MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
        MB.FullBasis{MB.Monomial,typeof(prod(x[1:2]))}(),
        tuple(),
    )
    @testset "Example 1 and 2" begin
        f = 1 + x[1]^4 * x[2]^2 + x[1]^2 * x[2]^4
        alg_el = MB.algebra_element(
            MB.sparse_coefficients(f),
            MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
        )
        with_var = SumOfSquares.Certificate.WithVariables(alg_el, x)
        newt = Certificate.NewtonDegreeBounds(tuple())
        @test SOS.Certificate.monomials_half_newton_polytope(
            monomials(f),
            newt,
        ) == [
            x[1]^2 * x[2],
            x[1] * x[2]^2,
            x[1]^2,
            x[1] * x[2],
            x[2]^2,
            x[1],
            x[2],
            1,
        ]
        @test SOS.Certificate.monomials_half_newton_polytope(
            monomials(f),
            Certificate.NewtonFilter(newt),
        ) == [x[1]^2 * x[2], x[1] * x[2]^2, 1]
        expected = Set(
            monomial_vector.([
                [x[1]^2 * x[2]],
                [x[1] * x[2]^2],
                [constant_monomial(x[1] * x[2])],
            ]),
        )
        for i in 0:2
            @test set_monos(
                Certificate.Sparsity.sparsity(
                    with_var,
                    Sparsity.Monomial(ChordalCompletion(), i),
                    certificate,
                ),
            ) == expected
        end
        expected = Set(
            monomial_vector.([
                [x[1]^2 * x[2]],
                [x[1] * x[2]^2, constant_monomial(x[1] * x[2])],
            ]),
        )
        @test set_monos(
            Certificate.Sparsity.sparsity(
                with_var,
                SignSymmetry(),
                certificate,
            ),
        ) == expected
    end
    @testset "Example 3 and 4" begin
        f = 1 + x[1]^4 + x[1] * x[2] + x[2]^4 + x[3]^2
        alg_el = MB.algebra_element(
            MB.sparse_coefficients(f),
            MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
        )
        with_var = SumOfSquares.Certificate.WithVariables(alg_el, x)
        @testset "$k $use_all_monomials" for k in 0:2,
            use_all_monomials in [false, true]

            if k == 1
                if use_all_monomials
                    @test set_monos(
                        Certificate.Sparsity.sparsity(
                            with_var,
                            Sparsity.Monomial(
                                ChordalCompletion(),
                                k,
                                use_all_monomials,
                            ),
                            certificate,
                        ),
                    ) == Set(
                        monomial_vector.([
                            [x[1]^2, x[2]^2, 1],
                            [x[1], x[2]],
                            [x[1] * x[2], 1],
                            [x[3]],
                        ]),
                    )
                else
                    @test set_monos(
                        Certificate.Sparsity.sparsity(
                            with_var,
                            Sparsity.Monomial(
                                ChordalCompletion(),
                                k,
                                use_all_monomials,
                            ),
                            certificate,
                        ),
                    ) == Set(
                        monomial_vector.([
                            [x[1]^2],
                            [x[2]^2],
                            [x[1], x[2]],
                            [x[1] * x[2], 1],
                            [x[3]],
                        ]),
                    )
                end
            else
                @test set_monos(
                    Certificate.Sparsity.sparsity(
                        with_var,
                        Sparsity.Monomial(
                            ChordalCompletion(),
                            k,
                            use_all_monomials,
                        ),
                        certificate,
                    ),
                ) == Set(
                    monomial_vector.([
                        [x[1], x[2]],
                        [x[3]],
                        [x[1]^2, x[2]^2, 1],
                        [x[1] * x[2], 1],
                    ]),
                )
            end
        end
        @test set_monos(
            Certificate.Sparsity.sparsity(
                SumOfSquares.Certificate.WithVariables(alg_el, x),
                SignSymmetry(),
                certificate,
            ),
        ) == Set(
            monomial_vector.([
                [x[1], x[2]],
                [x[3]],
                [x[1]^2, x[1] * x[2], x[2]^2, 1],
            ]),
        )
    end
end
function square_domain(ideal_certificate, d)
    @polyvar x y
    mult_cert = Certificate.MaxDegree(
        SOSCone(),
        MB.FullBasis{MB.Monomial,typeof(x * y)}(),
        MB.FullBasis{MB.Monomial,typeof(x * y)}(),
        d,
    )
    preorder_certificate =
        Certificate.Putinar(mult_cert, ideal_certificate(typeof(x * y)), d)
    f = x^2 * y^4 + x^4 * y^2 - 3 * x^2 * y * 2 + 1
    alg_el = MB.algebra_element(
        MB.sparse_coefficients(f),
        MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
    )
    K = @set(1 - x^2 >= 0 && 1 - y^2 >= 0)
    @testset "Square domain $k $use_all_monomials" for k in 0:4,
        use_all_monomials in [false, true]

        basis, preorder_bases = Certificate.Sparsity.sparsity(
            alg_el,
            K,
            Sparsity.Monomial(ChordalCompletion(), k, use_all_monomials),
            preorder_certificate,
        )
        if k == 1
            if use_all_monomials
                # [x, xy², x³], [1, y², x²],     [y, y³, x²y], [x, xy], [y, x², x²y], [1, x², x²y]
                # [x, xy², x³], [1, y, x², x²y], [y, y³, x²y], [x, xy], [y, x², x²y], [1, x², x²y]
                @test set_monos(basis) == Set(
                    monomial_vector.([
                        [x^2 * y, x^2, y],
                        [x^3, x * y^2, x],
                        [x^2 * y, x^2, 1],
                        [x^2, y^2, 1],
                        [x^2 * y, y^3, y],
                        [x * y, x],
                    ]),
                )
            else
                @test set_monos(basis) == Set(
                    monomial_vector.([
                        [x^2 * y, y^3],
                        [x^2, y],
                        [x^2 * y, 1],
                        [x^3, x * y^2],
                        [x * y, x],
                    ]),
                )
            end
        elseif k == 2
            @test set_monos(basis) == Set(
                monomial_vector.([
                    [x^2 * y, x^2, y^2, 1],
                    [x^3, x * y^2, x * y, x],
                    [x^2 * y, x^2, y, 1],
                    [x^2 * y, y^3, x^2, y],
                ]),
            )
        elseif k == 3
            @test set_monos(basis) == Set(
                monomial_vector.([
                    [x^3, x * y^2, x * y, x],
                    [x^2 * y, x^2, y^2, y, 1],
                    [x^2 * y, y^3, x^2, y, 1],
                ]),
            )
        else
            @test set_monos(basis) == Set(
                monomial_vector.([
                    [x^3, x * y^2, x * y, x],
                    [x^2 * y, y^3, x^2, y^2, y, 1],
                ]),
            )
        end
        expected = Set(monomial_vector.([[x^2, y^2, y, 1], [x * y, x]]))
        if k == 1
            if use_all_monomials
                @test set_monos(preorder_bases[1]) == Set(
                    monomial_vector.([[x * y, x], [x^2, y, 1], [x^2, y^2, 1]]),
                )
            else
                @test set_monos(preorder_bases[1]) == Set(
                    monomial_vector.([
                        [y, 1],
                        [x * y, x],
                        [x^2, y],
                        [x^2, y^2],
                    ]),
                )
            end
        else
            @test set_monos(preorder_bases[1]) == expected
        end
        if k == 1
            if use_all_monomials
                @test set_monos(preorder_bases[2]) == Set(
                    monomial_vector.([[x * y, x], [x^2, y], [x^2, y^2, 1]]),
                )
            else
                @test set_monos(preorder_bases[2]) == Set(
                    monomial_vector.([
                        [x * y, x],
                        [x^2, y],
                        [x^2, y^2],
                        [constant_monomial(x * y)],
                    ]),
                )
            end
        elseif k == 2
            @test set_monos(preorder_bases[2]) == Set(
                monomial_vector.([[x^2, y^2, 1], [x^2, y, 1], [x * y, x]]),
            )
        else
            @test set_monos(preorder_bases[2]) == expected
        end
    end
end
function sum_square(n)
    @testset "Sum square" begin
        @polyvar x[1:(2n)]
        certificate = Certificate.Newton(
            SOSCone(),
            MB.FullBasis{MB.Monomial,typeof(prod(x))}(),
            MB.FullBasis{MB.Monomial,typeof(prod(x))}(),
            tuple(),
        )
        f = sum((x[1:2:(2n-1)] .- x[2:2:(2n)]) .^ 2)
        alg_el = MB.algebra_element(
            MB.sparse_coefficients(f),
            MB.FullBasis{MB.Monomial,MP.monomial_type(f)}(),
        )
        expected = Set(
            monomial_vector.([
                monomial_vector([x[(2i-1)], x[2i], 1]) for i in 1:n
            ]),
        )
        @test set_monos(
            Certificate.Sparsity.sparsity(
                alg_el,
                Sparsity.Variable(),
                Certificate.MaxDegree(
                    SOSCone(),
                    MB.FullBasis{MB.Monomial,typeof(prod(x))}(),
                    MB.FullBasis{MB.Monomial,typeof(prod(x))}(),
                    2,
                ),
            ),
        ) == expected
        expected = Set(monomial_vector.([[x[(2i-1)], x[2i]] for i in 1:n]))
        @test set_monos(
            Certificate.Sparsity.sparsity(
                SumOfSquares.Certificate.WithVariables(alg_el, x),
                SignSymmetry(),
                certificate,
            ),
        ) == expected
    end
end
function drop_monomials()
    @testset "Drop monomials" begin
        @polyvar x
        basis = MB.SubBasis{MB.Monomial}([x^2])
        alg_el = _algebra_element([x^2])
        certificate = Certificate.MaxDegree(
            SOSCone(),
            MB.FullBasis{MB.Monomial,typeof(x^2)}(),
            MB.FullBasis{MB.Monomial,typeof(x^2)}(),
            2,
        )
        @testset "$k $use_all_monomials" for k in 0:2,
            use_all_monomials in [false, true]
            # The monomial `1˘ is dropped as it is useless.
            if use_all_monomials
                expected =
                    Set(monomial_vector.([[x], [constant_monomial(x^2)]]))
            else
                expected = Set([monomial_vector([x])])
            end
            @test set_monos(
                Certificate.Sparsity.sparsity(
                    SumOfSquares.Certificate.WithVariables(alg_el, [x]),
                    Sparsity.Monomial(
                        ChordalCompletion(),
                        k,
                        use_all_monomials,
                    ),
                    certificate,
                ),
            ) == expected
        end
        @testset "$(nameof(typeof(ideal_certificate)))" for ideal_certificate in
                                                            [
            Certificate.MaxDegree(
                SOSCone(),
                MB.FullBasis{MB.Monomial,typeof(x^2)}(),
                MB.FullBasis{MB.Monomial,typeof(x^2)}(),
                4,
            ),
            Certificate.Newton(
                SOSCone(),
                MB.FullBasis{MB.Monomial,typeof(x^2)}(),
                MB.FullBasis{MB.Monomial,typeof(x^2)}(),
                tuple(),
            ),
        ]
            preorder_certificate = Certificate.Putinar(
                Certificate.MaxDegree(
                    SOSCone(),
                    MB.FullBasis{MB.Monomial,typeof(x^2)}(),
                    MB.FullBasis{MB.Monomial,typeof(x^2)}(),
                    3,
                ),
                ideal_certificate,
                3,
            )
            K = @set x >= 0
            @testset "$k $use_all_monomials" for k in 0:3,
                use_all_monomials in [false, true]

                basis, preorder_bases = Certificate.Sparsity.sparsity(
                    _algebra_element([x^3]),
                    K,
                    Sparsity.Monomial(
                        ChordalCompletion(),
                        k,
                        use_all_monomials,
                    ),
                    preorder_certificate,
                )
                if ideal_certificate isa Certificate.Newton
                    if use_all_monomials
                        @test set_monos(basis) == Set(monomial_vector.([[x]]))
                        @test set_monos(preorder_bases[1]) ==
                              Set(monomial_vector.([[x, 1]]))
                    else
                        @test isempty(basis)
                        @test set_monos(preorder_bases[1]) ==
                              Set(monomial_vector.([[x]]))
                    end
                else
                    if k == 1 && !use_all_monomials
                        @test set_monos(basis) ==
                              Set(monomial_vector.([[x^2, x]]))
                    elseif (k == 2 && !use_all_monomials) ||
                           (k == 1 && use_all_monomials)
                        @test set_monos(basis) ==
                              Set(monomial_vector.([[x^2, 1], [x^2, x]]))
                    else
                        @test set_monos(basis) ==
                              Set(monomial_vector.([[x^2, x, 1]]))
                    end
                    if (k == 1 && !use_all_monomials)
                        @test set_monos(preorder_bases[1]) ==
                              Set(monomial_vector.([[x]]))
                    else
                        @test set_monos(preorder_bases[1]) ==
                              Set(monomial_vector.([[x, 1]]))
                    end
                end
            end
        end
    end
end
@testset "Sparsity" begin
    xor_complement_test()
    wml19()
    l09()
    square_domain(
        M -> Certificate.MaxDegree(
            SOSCone(),
            MB.FullBasis{MB.Monomial,M}(),
            MB.FullBasis{MB.Monomial,M}(),
            6,
        ),
        6,
    )
    square_domain(
        M -> Certificate.Newton(
            SOSCone(),
            MB.FullBasis{MB.Monomial,M}(),
            MB.FullBasis{MB.Monomial,M}(),
            tuple(),
        ),
        6,
    )
    sum_square(8)
    @test Certificate.Sparsity.appropriate_type(32) == Int64
    sum_square(32)
    @test Certificate.Sparsity.appropriate_type(64) == Int128
    sum_square(64)
    @test Certificate.Sparsity.appropriate_type(128) == BigInt
    sum_square(128)
    drop_monomials()
end
