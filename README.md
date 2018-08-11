# Sum of Squares Programming for Julia.

| **Documentation** | **PackageEvaluator** | **Build Status** | **Social** | **References to cite** |
|:-----------------:|:--------------------:|:----------------:|:----------:|:----------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] | [![Gitter][gitter-img]][gitter-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/en/a/af/Discourse_logo.png" width="64">][discourse-url] | |

**Important note**: The most recently tagged version of this package works with
most recently tagged version of [JuMP](https://github.com/JuliaOpt/JuMP.jl),
i.e. JuMP v0.18.x. The
in-development version of this package works with the in-development version of
[JuMP](https://github.com/JuliaOpt/JuMP.jl), i.e. JuMP 0.19-.

This packages contains the Sum of Squares reformulation for polynomial optimization.
When used in conjunction with [MultivariatePolynomial](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl), it provides a Sum of Squares Programming extension for JuMP.
Enabling the creation of sum of squares variables and constraints.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

Some presentations on, or using, SumOfSquares:
  * Benoît Legat at the JuMP Meetup 2017 [[Slides](http://www.juliaopt.org/meetings/mit2017/legat.pdf)] [[Video](https://youtu.be/kyo72yWYr54)]
  * [Joey Huchette at SIAM Opt 2017](https://docs.google.com/presentation/d/1ASfjB1TdLJmYxT0b6rnyGh9eLbMc-66bTOt3_3yvc90/edit?usp=sharing)

The following example shows how to find lower bounds for the Goldstein-Price function using this package with [MultivariatePolynomial](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl).

```julia
using MultivariatePolynomials
using JuMP
using PolyJuMP
using SumOfSquares
using DynamicPolynomials
using Mosek

# Create symbolic variables (not JuMP decision variables)
@polyvar x1 x2

# Create a Sum of Squares JuMP model with the Mosek solver
m = SOSModel(solver = MosekSolver())

# Create a JuMP decision variable for the lower bound
@variable m γ

# f(x) is the Goldstein-Price function
f1 = x1+x2+1
f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
f3 = 2*x1-3*x2
f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2

f = (1+f1^2*f2)*(30+f3^2*f4)

# Constraints f(x) - γ to be sum of squares
@constraint m f >= γ

@objective m Max γ

status = solve(m)

# The lower bound found is 3
println(getobjectivevalue(m))
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://juliaopt.github.io/SumOfSquares.jl/stable
[docs-latest-url]: https://juliaopt.github.io/SumOfSquares.jl/latest

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

[gitter-url]: https://gitter.im/JuliaOpt/SumOfSquares.jl?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/SumOfSquares.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt

[zenodo-url]: https://doi.org/10.5281/zenodo.1208672
[zenodo-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.1208672.svg
