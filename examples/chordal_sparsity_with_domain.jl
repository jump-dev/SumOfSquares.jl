using SumOfSquares
using MosekTools
using DynamicPolynomials

@polyvar x y z

K = @set 1-x^2-y^2>=0 && 1-y^2-z^2>=0

p = 1 - x*y - y*z

# use chordal sparsity

sparse_model = SOSModel(with_optimizer(Mosek.Optimizer))
@variable(sparse_model, t)
@objective(sparse_model, Max, t)

chordal_putinar(p - t, 2, K, model = sparse_model)

optimize!(sparse_model)
sparse_value = objective_value(sparse_model)

# no chordal sparsity

dense_model = SOSModel(with_optimizer(Mosek.Optimizer))
@variable(dense_model, t)
@objective(dense_model, Max, t)

@constraint(dense_model, p - t in SOSCone(), domain = K)

optimize!(dense_model)
dense_value = objective_value(dense_model)

println("Sparse and dense euqivalent?")
println("Absolute different of objective value: ", abs(dense_value - sparse_value))
