using SumOfSquares
using DynamicPolynomials
using MosekTools
using MathOptInterface
const MOI = MathOptInterface
using Profile, ProfileView
using Test

include("test_functions.jl")

function test(test_function, factory)
    m2 = SOSModel(factory)
    @variable m2 t
    @objective m2 Max t
    chordal_sos(test_function-t; model = m2)
    optimize!(m2)
    println(termination_status(m2))

    m = SOSModel(factory)
    @variable m t
    @objective m Max t
    @constraint m test_function-t in SOSCone()
    optimize!(m)
    println(termination_status(m))
    return abs(objective_value(m)-objective_value(m2))
end


@testset "choral_sos" begin
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
