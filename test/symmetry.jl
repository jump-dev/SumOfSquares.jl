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

function _test_orthogonal_transformation_to(A, B, As, Bs)
    U = SumOfSquares.Certificate.Symmetry.orthogonal_transformation_to(A, B, As, Bs)
    @test A ≈ U' * B * U
end

function _test_orthogonal_transformation_to(λ, As, Bs)
    A = sum(λ[i] * As[i] for i in eachindex(As))
    B = sum(λ[i] * Bs[i] for i in eachindex(Bs))
    _test_orthogonal_transformation_to(A, B, As, Bs)
    return
end

function _test_orthogonal_transformation_to(A, B)
    _test_orthogonal_transformation_to(A, B, typeof(A)[], typeof(B)[])
    _test_orthogonal_transformation_to(A, B, [A], [B])
    return
end

function _test_orthogonal_transformation_to(T::Type)
    A1 = T[
        0 -1
        1 0
    ]
    A2 = T[
        0 1
        -1 0
    ]
    D = Diagonal(T[1, -1])
    @test A1 == D * A2 * D
    _test_orthogonal_transformation_to(A1, A2)
    B1 = T[
        0 1
        1 0
    ]
    B2 = T[
        0 -1
        -1 0
    ]
    @test B1 == D * B2 * D
    _test_orthogonal_transformation_to(B1, B2)
    @test (A1 + B1) == D * (A2 + B2) * D
    _test_orthogonal_transformation_to(A1 + B1, A2 + B2)
    A1 = T[
        0 1
        -2 0
    ]
    A2 = T[
        0 -1
        2 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        0 -1
        -2 0
    ]
    A2 = T[
        0 1
        2 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        0 0 1
        1 0 0
        0 1 0
    ]
    A2 = T[
        0 1 0
        0 0 1
        1 0 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        1 0 0
        0 -1 0
        0 0 -1
    ]
    A2 = T[
        -1 0 0
        0 -1 0
        0 0 1
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        -1 0 1
        1 1 0
        0 1 -1
    ]
    A2 = T[
        -1 1 0
        0 1 1
        1 0 -1
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        -1 1 1
        2 1 0
        0 1 0
    ]
    A2 = T[
        0 1 0
        0 1 2
        1 1 -1
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[
        0 1 1
        2 0 0
        0 1 -1
    ]
    A2 = T[
        -1 1 0
        0 0 2
        1 1 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = ComplexF64[
        0 1 1
        2 0 0
        0 1 -1
    ]
    A2 = ComplexF64[
        -1 1 0
        0 0 2
        1 1 0
    ]
    _test_orthogonal_transformation_to(A1, A2)
    A1 = T[1 0; 0 -1]
    A2 = T[0 1; 1 0]
    _test_orthogonal_transformation_to([1, 1], [A1, A2], [-A1, A2])
    _test_orthogonal_transformation_to([2, 1], [A1, A2], [-A1, A2])
    _test_oeri_goluskin_orth(T)
    return
end

function _test_oeri_goluskin_orth(T)
    A1 = T[1 0; 0 1]
    A2 = T[-1 0; 0 -1]
    A3 = T[1 0; 0 -1]
    A4 = T[0.5 1
          1 0.5]
    A5 = T[-0.5 1
          1 -0.5]
    _test_orthogonal_transformation_to(ones(4), [A1, A2, A3, A4], [A1, A2, -A3, A5])
end

function test_orthogonal_transformation_to()
    @testset "$T" for T in [Int, Float64, ComplexF64]
        _test_orthogonal_transformation_to(T)
    end
end

function _test_block_diag(A, d)
    U = SumOfSquares.Certificate.Symmetry.ordered_block_diag(A, d)
    @test SumOfSquares.Certificate.Symmetry.ordered_block_check(U, A, d)
    return
end

function test_block_diag_dihedral()
    # From `dihedral.jl` example
    d = 2
    A2 = [
        0 1 0 0 0 0
        1 0 0 0 0 0
        0 0 0 0 0 1
        0 0 0 0 1 0
        0 0 0 1 0 0
        0 0 1 0 0 0
    ]
    A1 = [
        0 -1 0 0 0 0
        1 0 0 0 0 0
        0 0 0 0 0 -1
        0 0 0 0 1 0
        0 0 0 -1 0 0
        0 0 1 0 0 0
    ]
    _test_block_diag([A1, A2], d)
    # Using `GroupsCore.gens(G::DihedralGroup) = [DihedralElement(G.n, true, 1), DihedralElement(G.n, true, 0)]`, we get
    # see https://github.com/jump-dev/SumOfSquares.jl/issues/381#issuecomment-2296711306
    A1 = Matrix(Diagonal([1, -1, 1, -1, 1, -1]))
    _test_block_diag([A1, A2], d)
    return
end

function test_block_diag_alpha()
    α = 0.75
    A = [
        Matrix{Int}(I, 6, 6),
        [
            1 0 0 0 0 0
            0 -1 0 0 0 0
            0 0 0 0 1 0
            0 0 0 0 0 -1
            0 0 1 0 0 0
            0 0 0 -1 0 0
        ],
        [
            -0.5 α 0 0 0 0
            -α -0.5 0 0 0 0
            0 0 -0.5 α 0 0
            0 0 -α -0.5 0 0
            0 0 0 0 -0.5 α
            0 0 0 0 -α -0.5
        ],
    ]
    d = 2
    _test_block_diag(A, d)
    return
end

# See https://github.com/jump-dev/SumOfSquares.jl/issues/381
function test_oeri_goluskin()
    A1 = Matrix(1.0I, 16, 16)
    A2 = Matrix(-1.0I, 16, 16)
    A3 = Float64[
        1   0  0   0  0   0  0   0  0   0  0   0  0   0  0   0
        0  -1  0   0  0   0  0   0  0   0  0   0  0   0  0   0
        0   0  0   0  1   0  0   0  0   0  0   0  0   0  0   0
        0   0  0   0  0  -1  0   0  0   0  0   0  0   0  0   0
        0   0  1   0  0   0  0   0  0   0  0   0  0   0  0   0
        0   0  0  -1  0   0  0   0  0   0  0   0  0   0  0   0
        0   0  0   0  0   0  1   0  0   0  0   0  0   0  0   0
        0   0  0   0  0   0  0  -1  0   0  0   0  0   0  0   0
        0   0  0   0  0   0  0   0  1   0  0   0  0   0  0   0
        0   0  0   0  0   0  0   0  0  -1  0   0  0   0  0   0
        0   0  0   0  0   0  0   0  0   0  0   0  1   0  0   0
        0   0  0   0  0   0  0   0  0   0  0   0  0  -1  0   0
        0   0  0   0  0   0  0   0  0   0  1   0  0   0  0   0
        0   0  0   0  0   0  0   0  0   0  0  -1  0   0  0   0
        0   0  0   0  0   0  0   0  0   0  0   0  0   0  1   0
        0   0  0   0  0   0  0   0  0   0  0   0  0   0  0  -1
    ]
    A4 = sqrt.([
        1  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        3  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        0  0  0  0  1  3  0  0  0  0  0  0  0  0  0  0
        0  0  0  0  3  1  0  0  0  0  0  0  0  0  0  0
        0  0  1  3  0  0  0  0  0  0  0  0  0  0  0  0
        0  0  3  1  0  0  0  0  0  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  3  0  0  0  0  0  0  0  0
        0  0  0  0  0  0  3  1  0  0  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  1  3  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  3  1  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  0  0  0  0  1  3  0  0
        0  0  0  0  0  0  0  0  0  0  0  0  3  1  0  0
        0  0  0  0  0  0  0  0  0  0  1  3  0  0  0  0
        0  0  0  0  0  0  0  0  0  0  3  1  0  0  0  0
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  3
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  3  1
    ]) / 2
    d = 2
    _test_block_diag([A1, A2, A3, A4], d)
    return
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
