# SumOfSquares --- Sum of Squares Programming for Julia

[SumOfSquares](https://github.com/jump-dev/SumOfSquares.jl) implements Sum of Squares reformulation for [PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl),
enabling the creation of sum of squares variables and constraints in [JuMP](https://github.com/jump-dev/JuMP.jl).

The polynomial can be represented by any library implementing the [MultivariatePolynomial.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) interface.
That is, you can currently choose between [DynamicPolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [TypedPolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl).
As a rule of thumb, if you know at compile time (or at the time you are writing your code), the number of variable and that this number is small, use TypedPolynomials, otherwise, use DynamicPolynomials.

Some presentations on, or using, SumOfSquares (see [here](https://github.com/blegat/SumOfSquaresSlides) for the source code of the presentations):
  * Tillmann Weisser, Benoît Legat, Chris Coey, Lea Kapelevich and Juan Pablo Vielma at [JuliaCon 2019](https://juliacon.org/2019/) [[Slides](https://drive.google.com/open?id=1HiA-praFyejE0Z3nVSpFEv938TAcPjA9)] [[Video](https://www.youtube.com/watch?v=cTmqmPcroFo)]
  * Benoît Legat at [CNLS 2019](https://cnls.lanl.gov/External/showtalksummary.php?selection=7768) [[Slides](https://drive.google.com/open?id=1kNF18C7RY2zi7jcZBMO1PRXtHuvVTFPn)]
  * Benoît Legat at [EURO 2019](https://www.euro2019dublin.com/) [[Slides](https://drive.google.com/open?id=1Wry56NzzL4QBRSwuhP4AlKOe2i2FL7dk)]
  * Benoît Legat at [juliaday Nantes 2019](https://julialang.univ-nantes.fr/programme/) [[Slides](https://drive.google.com/open?id=1pN3G9Pr8jbzK9EEaJ9a6p_qKwSbxb2bo)]
  * Benoît Legat at [Summer School on Numerical Computing in Algebraic Geometry 2018](https://www.mis.mpg.de/calendar/conferences/2018/nc2018.html) [[Poster](https://drive.google.com/open?id=1pf9rdoVEjAnD164rptLki1AG0AH4i88M)]
  * Benoît Legat at [The First Annual JuMP-dev Workshop 2017](https://jump.dev/meetings/mit2017/) [[Slides](https://drive.google.com/file/d/1ea5eSMvMB3jXPuljzNGmMKied-n50YIo/view?usp=sharing)] [[Video](https://youtu.be/kyo72yWYr54)]
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
