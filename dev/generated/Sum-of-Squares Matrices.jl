using Test
using DynamicPolynomials
@polyvar x
P = [x^2 - 2x + 2 x
            x     x^2]

using SumOfSquares

using CSDP
factory = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = SOSModel(factory)
mat_cref = @constraint(model, P in PSDCone())
optimize!(model)
@test termination_status(model) == MOI.OPTIMAL

@test length(certificate_monomials(mat_cref)) == 3

@polyvar y[1:2]
p = vec(y)' * P * vec(y)
nothing #md See https://github.com/JuliaDocs/Documenter.jl/issues/1387

X = monomials(p)
@test Certificate.monomials_half_newton_polytope(X, tuple(), apply_post_filter = false) == [x * y[1], x * y[2], y[1] * y[2], x, y[1], y[2]]

@test Certificate.monomials_half_newton_polytope(X, ([x], y), apply_post_filter = false) == [x * y[1], x * y[2], y[1], y[2]]

@test Certificate.monomials_half_newton_polytope(X, ([x], y)) == [x * y[1], x * y[2], y[1]]

@test Certificate.monomials_half_newton_polytope(X, tuple()) == [x * y[1], x * y[2], y[1]]

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

