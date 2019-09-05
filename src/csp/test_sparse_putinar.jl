using SumOfSquares
using DynamicPolynomials
using MosekTools
using Profile, ProfileView
using Test

include("src/csp/test_functions.jl")

@testset "chained_singular" begin
    n = 100
    m = SOSModel(with_optimizer(Mosek.Optimizer))
    @variable m t
    @objective m Max t
    chordal_sos(chained_singular(n)-t; model = m)
    optimize!(m)
end

function test(n)
    m = SOSModel(with_optimizer(Mosek.Optimizer))
    @variable m t
    @objective m Max t
    chordal_sos(chained_singular(n)-t; model = m)
    optimize!(m)
end


@testset "broyden_banded" begin
end

@testset "broyden_tridiagonal" begin
end

@testset "chained_wood" begin
end

@testset "generalized_rosenbrock" begin
end


