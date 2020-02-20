using SumOfSquares
using MosekTools
using DynamicPolynomials


function sos_lower_bound(factory, sparse::Sparsity)
    @polyvar x y z
    K = @set 1-x^2-y^2>=0 && 1-y^2-z^2>=0
    p = 1 - x*y - y*z
    model = Model(factory)
    @variable(model, t)
    @objective(model, Max, t)
    @constraint(model, p - t in SOSCone(), domain = K, sparse = sparse)
    optimize!(model)
    return objective_value(model)
end

const factory = Mosek.Optimizer

# use chordal sparsity
sparse_value = sos_lower_bound(factory, VariableSparsity())

# no chordal sparsity
dense_value = sos_lower_bound(factory, NoSparsity())

println("Sparse and dense equivalent?")
println("Absolute different of objective value: ", abs(dense_value - sparse_value))
