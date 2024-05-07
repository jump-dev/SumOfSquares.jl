using SumOfSquares
using DynamicPolynomials

@polyvar x[1:2]

import Clarabel
using Dualization
model = SOSModel(dual_optimizer(Clarabel.Optimizer))

@variable(model, γ)

f1 = x[1] + x[2] + 1
f2 = 19 - 14*x[1] + 3*x[1]^2 - 14*x[2] + 6*x[1]*x[2] + 3*x[2]^2
f3 = 2*x[1] - 3*x[2]
f4 = 18 - 32*x[1] + 12*x[1]^2 + 48*x[2] - 36*x[1]*x[2] + 27*x[2]^2
f = (1 + f1^2*f2) * (30 + f3^2*f4)

con_ref = @constraint(model, f >= γ)
@objective(model, Max, γ)
optimize!(model)

solution_summary(model)

ν = moment_matrix(con_ref)

μ = measure(ν, atol = 1e-5)

ν_truncated = moment_matrix(μ, monomials(x, 0:3))

using LinearAlgebra
LinearAlgebra.svdvals(Matrix(ν_truncated.Q))
LinearAlgebra.rank(Matrix(ν_truncated.Q), rtol = 1e-3)
svdvals(Matrix(ν_truncated.Q))

svdvals(Matrix(ν.Q))

atomic_measure(ν, FixedRank(3))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
