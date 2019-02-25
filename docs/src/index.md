# SumOfSquares --- Sum of Squares Programming for Julia

[SumOfSquares](https://github.com/JuliaOpt/SumOfSquares.jl) implements Sum of Squares reformulation for [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl),
enabling the creation of sum of squares variables and constraints in [JuMP](https://github.com/JuliaOpt/JuMP.jl).

The polynomial can be represented by any library implementing the [MultivariatePolynomial.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) interface.
That is, you can currently choose between [DynamicPolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [TypedPolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl).
As a rule of thumb, if you know at compile time (or at the time you are writing your code), the number of variable and that this number is small, use TypedPolynomials, otherwise, use DynamicPolynomials.

Some presentations on, or using, SumOfSquares:
  * Benoît Legat at the JuMP Meetup 2017 [[Slides](http://www.juliaopt.org/meetings/mit2017/legat.pdf)] [[Video](https://youtu.be/kyo72yWYr54)]
  * [Joey Huchette at SIAM Opt 2017](https://docs.google.com/presentation/d/1ASfjB1TdLJmYxT0b6rnyGh9eLbMc-66bTOt3_3yvc90/edit?usp=sharing)

The following example shows how to find lower bounds for the Goldstein-Price function using this package.

```julia
using SumOfSquares
using DynamicPolynomials
using MosekTools

# Create symbolic variables (not JuMP decision variables)
@polyvar x1 x2

# Create a Sum of Squares JuMP model with the Mosek solver
model = SOSModel(with_optimizer(Mosek.Optimizer))

# Create a JuMP decision variable for the lower bound
@variable(model, γ)

# f(x) is the Goldstein-Price function
f1 = x1+x2+1
f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
f3 = 2*x1-3*x2
f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2

f = (1+f1^2*f2)*(30+f3^2*f4)

# Constraints f(x) - γ to be sum of squares
@constraint(model, f >= γ)

@objective(model, Max, γ)

optimize!(model)

# The lower bound found is 3
println(objective_value(model))
```

## Contents
```@contents
Pages = ["sumofsquares.md", "variables.md", "constraints.md"]
Depth = 2
```
