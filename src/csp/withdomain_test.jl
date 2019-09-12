using SumOfSquares
using MosekTools
using DynamicPolynomials

@polyvar x y z

K = @set 1-x^2-y^2>=0 && 1-y^2-z^2>=0

p = 1 - x*y - y*z

# use chordal sparsity

sm = SOSModel(with_optimizer(Mosek.Optimizer))
@variable sm t
@objective sm Max t
chordal_putinar(p-t, 2, K,  model = sm)

optimize!(sm)
sparse_value = objective_value(sm)


# no chordal sparsity

m = SOSModel(with_optimizer(Mosek.Optimizer))
@variable m t
@objective m Max t

@constraint m p-t in SOSCone() domain = K 

optimize!(m)
dense_value = objective_value(m)

println("Sparse and dense euqivalent?")
abs(dense_value-sparse_value)<10^-8
