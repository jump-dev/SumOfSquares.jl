module TestSymmetry

using LinearAlgebra, SparseArrays, Test
using SumOfSquares

function test_linsolve()
    x = [1, 2]
    for A in [
        [1 0 2 3
         0 1 3 -2],
        [1 2 0 3
         0 3 1 -2],
    ]
        b = A' * x
        @test Certificate.Symmetry._linsolve(A, b, Symmetry._RowEchelonMatrix()) ≈ x
        B = float.(A)
        for i in 1:size(B, 1)
            B[i, :] = normalize(B[i, :])
        end
        b = B' * x
        @test Certificate.Symmetry.__linsolve(B, b, Symmetry._OrthogonalMatrix()) ≈ x
        @test Certificate.Symmetry.__linsolve(sparse(B), b, Symmetry._OrthogonalMatrix()) ≈ x
    end
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
