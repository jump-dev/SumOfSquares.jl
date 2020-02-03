using Test
using SumOfSquares
using DynamicPolynomials

function rearrangement_test(optimizer, config::MOIT.TestConfig)
    atol = config.atol
    rtol = config.rtol

    model = _model(optimizer)

    @polyvar x y z
    # If x, y is in the same order as y, z, then by the rearrangement inequality
    # xy + yz ≤ xz + y^2
    # Then, if x and z are nonnegative, by the Geometric Mean-Quadratic Mean inequality,
    # xz ≤ (x^2 + z^2)/2
    # Therefore,
    # 2xy + 2yz ≤ x^2 + z^2 + 2y^2 over the set:
    S = @set y ≤ x && y ≤ z && x ≥ 0 && z ≥ 0
    # In fact, it is nonnegative everywhere so let's add terms `x + z`:
    cref = @constraint(model, x^2 + z^2 + 2y^2 - 2x * y - 2y * z + x + z in SOSCone(),
                       domain = S, sparse = VariableSparsity())

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL

    @test primal_status(model) == MOI.FEASIBLE_POINT

    p = gram_matrix(cref)
    @test p isa SumOfSquares.SparseGramMatrix
    @test length(p.sub_gram_matrices) == 2
    @test getmat(p.sub_gram_matrices[1]) ≈ [1 -1 0; -1 1 0; 0 0 0] atol=9atol rtol=9rtol
    @test p.sub_gram_matrices[1].basis.monomials == [y, z, 1]
    @test getmat(p.sub_gram_matrices[2]) ≈ [1 -1 0; -1 1 0; 0 0 0] atol=9atol rtol=9rtol
    @test p.sub_gram_matrices[2].basis.monomials == [x, y, 1]

    λ = lagrangian_multipliers(cref)
    @test length(λ) == 4
    @test λ[1] isa SumOfSquares.SparseGramMatrix
    @test length(λ[1].sub_gram_matrices) == 1
    @test λ[1].sub_gram_matrices[1].basis.monomials == [x, y, 1]
    @test λ[2] isa SumOfSquares.SparseGramMatrix
    @test length(λ[2].sub_gram_matrices) == 1
    @test λ[2].sub_gram_matrices[1].basis.monomials == [y, z, 1]
    @test λ[3] isa SumOfSquares.SparseGramMatrix
    @test length(λ[3].sub_gram_matrices) == 1
    @test λ[3].sub_gram_matrices[1].basis.monomials == [x, y, 1]
    @test λ[4] isa SumOfSquares.SparseGramMatrix
    @test length(λ[4].sub_gram_matrices) == 1
    @test λ[4].sub_gram_matrices[1].basis.monomials == [y, z, 1]

    ν = moment_matrix(cref)
    @test ν isa SparseMomentMatrix
    @test length(ν.sub_moment_matrices) == 2
    @test ν.sub_moment_matrices[1].basis.monomials == [y, z, 1]
    @test ν.sub_moment_matrices[2].basis.monomials == [x, y, 1]
end
sd_tests["rearrangement"] = rearrangement_test
