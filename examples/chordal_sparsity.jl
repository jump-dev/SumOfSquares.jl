using SumOfSquares
using DynamicPolynomials
using Test

include("sparse_polynomials.jl")

function test(test_function, optimizer_constructor)
    sparse_model = sos_lower_bound(test_function, optimizer_constructor, Sparsity.Variable())
    dense_model  = sos_lower_bound(test_function, optimizer_constructor, Sparsity.NoPattern())
    return abs(objective_value(sparse_model) - objective_value(dense_model))
end

function chordal_sos(optimizer_constructor, tol=1e-5)
    val = test(chained_singular(4), optimizer_constructor)
    @test val < tol

    val = test(broyden_banded(2), optimizer_constructor)
    @test val < tol

    val = test(broyden_tridiagonal(4), optimizer_constructor)
    @test val < tol

    val = test(chained_wood(4), optimizer_constructor)
    @test val < tol

    val = test(generalized_rosenbrock(4), optimizer_constructor)
    @test val < tol
end

import CSDP
chordal_sos(CSDP.Optimizer)
