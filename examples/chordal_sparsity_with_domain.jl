using SumOfSquares
using DynamicPolynomials
using Test

function sos_lower_bound(factory, sparsity::Sparsity)
    @polyvar x y z
    K = @set 1 - x^2 - y^2 >= 0 && 1 - y^2 - z^2 >= 0
    p = 1 - x*y - y*z
    model = Model(factory)
    @variable(model, t)
    @objective(model, Max, t)
    @constraint(model, p - t in SOSCone(), domain = K, sparsity = sparsity)
    optimize!(model)
    return objective_value(model)
end

function sparse_vs_dense(optimizer_constructor)
    # use chordal sparsity
    sparse_value = sos_lower_bound(optimizer_constructor, VariableSparsity())

    # no chordal sparsity
    dense_value = sos_lower_bound(optimizer_constructor, NoSparsity())

    @test abs(sparse_value - dense_value) < 1e-6
end

import CSDP
sparse_vs_dense(CSDP.Optimizer)
