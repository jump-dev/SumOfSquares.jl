# Polynomial and Sum-of-Squares variables

## Polynomial variables

While [JuMP](https://github.com/JuliaOpt/JuMP.jl) allows to create decision
variables representing a number whose value needs to be optimized upon by the
optimizer, [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) allows to create
*polynomial* decision variables. In order to do that, we first need to create
polynomial variables with the `@polyvar` macro:
```julia
using DynamicPolynomials # or TypedPolynomials, pick your favorite
@polyvar x y
```
Note that these should not be confused with JuMP's decision variables which are
created using the `@variable` macro. The polynomial decision variable that will
be created need to be parametrized by a polynomial basis of finite size.
For instance, if we want to find a quadratic polynomial, we can parametrize it
with all monomials of degree between 0 and 2. Generating a vector with such
monomials can be achieved through the `monomials` function:
```julia
using MultivariatePolynomials
X = monomials([x, y], 0:2)
```
We can now create our polynomial variable `p` as follows:
```julia
using PolyJuMP
@variable(model, p, Poly(X))
```
This creates a vector of decision variables `a` and sets `p` as the scalar
product between `a` and `X`.

Just like with classical JuMP's decision variables, containers of polynomial
variables can be created as follows:
```julia
using PolyJuMP
@variable(model, p[1:3, 1:4], Poly(X))       # Creates a Matrix
@variable(model, p[[:a, :b], -2:2], Poly(X)) # Creates a JuMPArray
@variable(model, p[i=1:3, j=i:3], Poly(X))   # Creates a Dict
```

For more flexibility, polynomials parametrized by decision variables can also
be created "by hand" for instance as follows:
```julia
@variable(model, α)
@variable(model, β)
p = α*x^2 + (α+β)*y^2*x + β*y^3
```
The coefficients of the polynomial can even be quadratic, e.g.
```julia
@variable(model, α)
@variable(model, β)
p = (3α^2+β)*x^2 + (α*β+2β)*y^2*x + β*y^3
```
Let me stress again the distinction between `α` and `β` which are *decision*
variables and `x` and `y` which are *polynomial* variables.

## Nonnegative polynomial variables

In order to create a sum-of-squares polynomial variable, the syntax is exactly
the same except `SOSPoly` should be used instead of `Poly`.
For instance, the following code creates a ``3 \\times 4`` matrix of
sum-of-squares polynomial variables:
```julia
using SumOfSquares
@variable(model, p[1:3, 1:4], SOSPoly(X))
```
There is however an *important* difference between the signification of the
vector of monomials `X` between `Poly` and `SOSPoly`. For `SOSPoly`, it
creates a positive semidefinite matrix of variables `Q` and sets `p` as the
value of `X' * Q * X`. That is, for instance, if `X` contains all the monomials
of degree 2, then all the monomials of `p` will have degree 4 (i.e. `p` will be
a quartic form).

Similarly, to create diagonally-dominant-sum-of-squares polynomial variables
(see [Definition 2, AM17]), use `DSOSPoly(X)`. This creates a diagonally
dominant matrix of variables `Q` and sets the polynomial variables as the value
of `X' * Q * X`.

Finally, to create scaled-diagonally-dominant-sum-of-squares polynomial
variables (see [Definition 2, AM17]), use `DSOSPoly(X)`. This creates a
scaled diagonally dominant matrix of variables `Q` and sets the polynomial
variables as the value of `X' * Q * X`.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, 2017
