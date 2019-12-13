using SumOfSquares
using DynamicPolynomials
using MosekTools
using Test

include("sparse_polynomials.jl")

function test(test_function, factory)
    sparse_model = sparse(test_function, factory)
    dense_model = dense(test_function, factory)
    return abs(objective_value(sparse_model) - objective_value(dense_model))
end

@testset "chordal_sos" begin
    factory = with_optimizer(Mosek.Optimizer)

    val = test(chained_singular(12), factory)
    @test abs(val) < 1e-6

    val = test(broyden_banded(6), factory)
    @test abs(val) < 1e-6

    val = test(broyden_tridiagonal(12), factory)
    @test abs(val) < 1e-6

    val = test(chained_wood(12), factory)
    @test abs(val) < 1e-6

    val = test(generalized_rosenbrock(12), factory)
    @test abs(val) < 1e-6

end
