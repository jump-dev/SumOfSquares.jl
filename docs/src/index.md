# SumOfSquares --- Sum of Squares Programming for Julia

[SumOfSquares](https://github.com/jump-dev/SumOfSquares.jl) implements Sum of Squares reformulation for [PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl),
enabling the creation of sum of squares variables and constraints in [JuMP](https://github.com/jump-dev/JuMP.jl).

The polynomial can be represented by any library implementing the [MultivariatePolynomial.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) interface.
That is, you can currently choose between [DynamicPolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [TypedPolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl).
As a rule of thumb, if you know at compile time (or at the time you are writing your code), the number of variable and that this number is small, use TypedPolynomials, otherwise, use DynamicPolynomials.

## Contents
```@contents
Pages = ["sumofsquares.md", "variables.md", "constraints.md"]
Depth = 2
```
