
using SumOfSquares
using MultivariatePolynomials
const MP = MultivariatePolynomials
using DynamicPolynomials
const CEG = SumOfSquares.Certificate.CEG
using MosekTools
include("../src/Certificate/MonoSparse/monosparse.jl")
using MosekTools

d = 6
@polyvar x y 
f = x^2*y^4 + x^4*y^2 - 3*x^2*y*2 + 1
K = @set(1-x^2>=0 && 1-y^2>=0)

println("Dense Model:")
m = SOSModel(with_optimizer(Mosek.Optimizer))
@variable m t
@objective m Max t
@constraint m f-t >= 0 domain = K
optimize!(m)
println(termination_status(m))

println()
println("MonoSparse Model:")
sm = SOSModel(with_optimizer(Mosek.Optimizer))
@variable sm st
@objective sm Max st
mono_sparse_certificate(sm, f-st, K, d)
optimize!(sm)
println(termination_status(sm))

println("Difference objective values:")
println(objective_value(m) - objective_value(sm))
