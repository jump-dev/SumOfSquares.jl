module TestSymmetry

using LinearAlgebra, SparseArrays, Test
using SumOfSquares

function test_linsolve()
    x = [1, 2]
    for A in [
        [
            1 0 2 3
            0 1 3 -2
        ],
        [
            1 2 0 3
            0 3 1 -2
        ],
    ]
        b = A' * x
        @test Certificate.Symmetry._linsolve(
            A,
            b,
            Symmetry._RowEchelonMatrix(),
        ) ≈ x
        B = float.(A)
        for i in axes(B, 1)
            B[i, :] = normalize(B[i, :])
        end
        b = B' * x
        @test Certificate.Symmetry.__linsolve(
            B,
            b,
            Symmetry._OrthogonalMatrix(),
        ) ≈ x
        @test Certificate.Symmetry.__linsolve(
            sparse(B),
            b,
            Symmetry._OrthogonalMatrix(),
        ) ≈ x
    end
end

function _test_orthogonal_transformation_to(A, B)
    U = SumOfSquares.Certificate.Symmetry.orthogonal_transformation_to(A, B)
    @test A ≈ U' * B * U
end

function test_orthogonal_transformation_to()
    A1 = [
        0 -1
        1 0
    ]
    A2 = [
        0 1
        -1 0
    ]
    D = Diagonal([1, -1])
    @test A1 ≈ D * A2 * D
    _test_orthogonal_transformation_to(A1, A2)
    B1 = [
        0 1
        1 0
    ]
    B2 = [
        0 -1
        -1 0
    ]
    @test B1 ≈ D * B2 * D
    _test_orthogonal_transformation_to(B1, B2)
    @test (A1 + B1) ≈ D * (A2 + B2) * D
    _test_orthogonal_transformation_to(A1 + B1, A2 + B2)
    A1 = [
        0 1
        -2 0
    ]
    A2 = [
        0 -1
        2 0
    ]
    return _test_orthogonal_transformation_to(A1, A2)
end

function test_block_diag()
    # From `dihedral.jl` example
    A = [
        [
            0 -1 0 0 0 0
            1 0 0 0 0 0
            0 0 0 0 0 -1
            0 0 0 0 1 0
            0 0 0 -1 0 0
            0 0 1 0 0 0
        ],
        [
            0 1 0 0 0 0
            1 0 0 0 0 0
            0 0 0 0 0 1
            0 0 0 0 1 0
            0 0 0 1 0 0
            0 0 1 0 0 0
        ],
    ]
    d = 2
    U = SumOfSquares.Certificate.Symmetry.ordered_block_diag(A, d)
    @test SumOfSquares.Certificate.Symmetry.ordered_block_check(U, A, d)
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestSymmetry.runtests()
