using DynamicPolynomials
@polyvar x y
motzkin = x^4*y^2 + x^2*y^4 + 1 - 3x^2*y^2

using SumOfSquares
import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = SOSModel(solver)
@constraint(model, motzkin >= 0) # We constraint `motzkin` to be a sum of squares

optimize!(model)

termination_status(model)

dual_status(model)

model = SOSModel(solver)
@constraint(model, (x^2 + y^2) * motzkin >= 0) # We constraint the `(x^2 + y^2) * motzkin` to be a sum of squares

optimize!(model)

termination_status(model)

primal_status(model)

model = SOSModel(solver)
X = monomials([x, y], 0:2)

@variable(model, deno, Poly(X))

@constraint(model, deno >= 1)
@constraint(model, deno * motzkin >= 0)
optimize!(model)

termination_status(model)

primal_status(model)

value(deno)

using RecipesBase
@recipe function f(x::AbstractVector, y::AbstractVector, p::Polynomial)
    x, y, (x, y) -> p(variables(p) => [x, y])
end
import Plots
Plots.plot(
    range(-2, stop=2, length=100),
    range(-2, stop=2, length=100),
    motzkin,
    st = [:surface],
    seriescolor=:heat,
    colorbar=:none,
    clims = (-10, 80)
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
