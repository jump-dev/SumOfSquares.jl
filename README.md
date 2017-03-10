# Sum of Squares Programming for Julia.

| **PackageEvaluator** | **Build Status** |
|:--------------------:|:----------------:|
| [![][pkg-0.5-img]][pkg-0.5-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![][pkg-0.6-img]][pkg-0.6-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

This packages contains the Sum of Squares reformulation for polynomial optimization.
When used in conjunction with [MultivariatePolynomial.jl](https://github.com/blegat/MultivariatePolynomials.jl) and [PolyJuMP.jl](https://github.com/JuliaOpt/PolyJuMP.jl), it provides a Sum of Squares Programming extension for JuMP.
Enabling the creation of sum of squares variables and constraints.

The following example shows how to find lower bounds for the Goldstein-Price function using this package with [MultivariatePolynomial.jl](https://github.com/blegat/MultivariatePolynomials.jl) and [PolyJuMP.jl](https://github.com/JuliaOpt/PolyJuMP.jl).

```julia
using MultivariatePolynomials
using JuMP
using PolyJuMP
using SumOfSquares

# Create symbolic variables (not JuMP decision variables)
@polyvar x1 x2

# Create a JuMP model with the default SDP solver (you should have at least one installed)
m = Model()

# Create a JuMP decision variable for the lower bound
@variable m γ

# f(x) is the Goldstein-Price function
f1 = x1+x2+1
f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
f3 = 2*x1-3*x2
f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2

f = (1+f1^2*f2)*(30+f3^2*f4)

# Constraints f(x) - γ to be sum of squares
@polyconstraint m f >= γ

@objective m Max γ

status = solve(m)

# The lower bound found is 3
println(getobjectivevalue(m))
```

[pkg-0.5-img]: http://pkg.julialang.org/badges/SumOfSquares_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=SumOfSquares
[pkg-0.6-img]: http://pkg.julialang.org/badges/SumOfSquares_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=SumOfSquares

[build-img]: https://travis-ci.org/JuliaOpt/SumOfSquares.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaOpt/SumOfSquares.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/o49y96hl1xl5aytn?svg=true
[winbuild-url]: https://ci.appveyor.com/project/JuliaOpt/sumofsquares-jl
[coveralls-img]: https://coveralls.io/repos/github/JuliaOpt/SumOfSquares.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/JuliaOpt/SumOfSquares.jl?branch=master
[codecov-img]: http://codecov.io/github/JuliaOpt/SumOfSquares.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaOpt/SumOfSquares.jl?branch=master
