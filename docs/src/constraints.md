# Constraints

## Equality constraints between polynomials

Equality between polynomials in
[PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) uses the same syntax as
equality between affine or quadratic expression in
[JuMP](https://github.com/JuliaOpt/JuMP.jl).
For instance, creating two quadratic `n`-variate polynomials `p` and `q` that
must sum up to one can be done as follows:
```julia
using DynamicPolynomials
@polyvar x[1:n]
using MultivariatePolynomials
X = monomials(x, 0:2)
using PolyJuMP
@variable(model, p, Poly(X))
@variable(model, q, Poly(X))
@constraint(model, p + q == 1)
```

Vectorized constraints can also be used as well as vector of constraints,
named constraints, ...
For instance instance, if `P` and `Q` are two ``n \times n`` matrices of
polynomials, the following constraints the sum of rows and columns to match
```julia
@constraint(model, con[i=1:n], sum(P[i, :]) == sum(Q[:, i]))
```
and `con[i]` contains the reference to the constraint between the `i`th row
of `P` and the `i`th column of `Q`.

## Inequality constraints between polynomials

Polynomials can be constrained to be sum-of-squares with the `in` syntax.
For instance, to constraints a polynomial `p` to be sum-of-squares, do
```julia
@constraint(model, p in SOSCone())
```

### Automatically interpreting polynomial nonnegativity as a sum-of-squares constraint

As detailed in [When is nonnegativity equivalent to sum of squares ?](@ref),
the nonnegativity of a polynomial is not equivalent to the existence of a
sum-of-squares decomposition. However, if explicitely specified, nonnegativity
constraints can be automatically interpreted as sum-of-squares constraints.
The simplest way to do that is to create the model with
```julia
model = SOSModel(...)
```
instead of
```julia
model = Model(...)
```
An alternative equivalent way is to call `setpolymodule!` after creating the
model:
```julia
setpolymodule!(model, SumOfSquares)
```
This second approach may be useful if the SumOfSquares JuMP extension need to
be used with another JuMP extension that also has a special model constructor.
A third alternative is the following:
```julia
PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, SOSCone)
PolyJuMP.setdefault!(model, PolyJuMP.NonNegPolyMatrix, SOSMatrixCone)
```
This approach adds the flexibility to choose the default cone for

* constraints of the form
  `@constraint(mode, ..., some_polynomial ≥ other_polynomial, ...)`
  which is the cone given as default to `PolyJuMP.NonNegPoly`; and
* constraints of the form
  `@constraint(mode, ..., some_matrix_of_polynomial in PSDCone(), ...)`
  or
  `@SDconstraint(mode, ..., some_matrix_of_polynomial ⪰ other_matrix_of_polynomial, ...)`
  which is the cone given as default to `PolyJuMP.NonNegPolyMatrix`.

For instance, to use the diagonally-dominant-sum-of-squares cone (see
[Definition 2, AM17]) for the first type of contraints, do
```julia
PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, DSOSCone)
```
## Changing the polynomial basis

As introduced in [Choosing a polynomial basis](@ref), there may be numerical
advantages to use another basis than the standard monomial basis when creating
polynomial variables. Similarly, other polynomial basis can be used for
polynomial constraints. However, for constraints, the polynomial space is
determined by the polynomial constrained to be nonnegative. For instance,
consider the constraint:
```julia
@constraint(model, α * x^2 + β * y^2 ≥ (α - β) * x * y)
```
where `α` and `β` are JuMP decision variables and `x` and `y` are polynomial
variables. Since the polynomial is a quadratic form, the sum-of-squares
certificate is also a quadratic form (see [Section~3.3.4, BPT12]). Hence the
default polynomial basis used for the [Nonnegative polynomial variables]
certificate is `MonomialBasis([x, y])`, that is, we search for a positive
semidefinite matrix `Q` such that
```math
α x^2 + β y^2 - (α - β) x y = X^\top Q X
```
where ``X = (x, y)``.

As the polynomial space is determined by the polynomial being constrained,
only the basis *type* need to be given.For instance, to use the scaled monomial
basis, use
```julia
@constraint(model, p ≥ q, basis = ScaledMonomialBasis)
```

## Polynomial nonnegativity on a subset of the space

By default, the constraint
```julia
@constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α)
```
constrains the polynomial to be nonnegative for every real numbers `x` and `y`.
However, the set of points `(x, y)` for which the polynomial is constrained
to be nonnegative can be specified by the `domain` keyword:
```julia
using SemialgebraicSets
S = @set x >= 0 && y >= 0 && x + y >= 1
@constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α, domain = S)
```
See [this notebook](https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb)
for a detailed example.

## Dual of polynomial constraints

The dual of a polynomial constraint `cref` is a moment serie `μ` as defined in
[MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl).
The dual can be obtained with the `JuMP.resultdual` function as with classical
dual values in JuMP. The matrix of moment can be obtaine as follows:
```julia
μ = JuMP.resultdual(cref)
ν = matmeasure(μ, certificate_monomials(cref))
```
The `extractatoms` function of [MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl)
can be used to check if there exists an *atomic* measure (i.e. a measure that is
a sum of Dirac measures) for that has the moments given in `ν`.
This can be used for instance in polynomial optimization (see
[this notebook](https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb))
or stability analysis (see
[this notebook](https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/LPJ17e43.ipynb)).

### References

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, 2012.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, 2017
