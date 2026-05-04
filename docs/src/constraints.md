# Constraints

## Equality constraints between polynomials

Equality between polynomials in
[PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl) uses the same syntax as
equality between affine or quadratic expression in
[JuMP](https://github.com/jump-dev/JuMP.jl).
For instance, creating two quadratic `n`-variate polynomials `p` and `q` that
must sum up to one can be done as follows:
```jldoctest constraint-pq
julia> n = 3
3

julia> using DynamicPolynomials

julia> @polyvar x[1:n]
(Variable{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}[x₁, x₂, x₃],)

julia> X = monomials(x, 0:2)
10-element MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}:
 1
 x₃
 x₂
 x₁
 x₃²
 x₂x₃
 x₂²
 x₁x₃
 x₁x₂
 x₁²

julia> using SumOfSquares

julia> model = Model();

julia> @variable(model, p, Poly(X))
(_[1])·1 + (_[2])·x₃ + (_[3])·x₂ + (_[4])·x₁ + (_[5])·x₃² + (_[6])·x₂x₃ + (_[7])·x₂² + (_[8])·x₁x₃ + (_[9])·x₁x₂ + (_[10])·x₁²

julia> @variable(model, q, Poly(X))
(_[11])·1 + (_[12])·x₃ + (_[13])·x₂ + (_[14])·x₁ + (_[15])·x₃² + (_[16])·x₂x₃ + (_[17])·x₂² + (_[18])·x₁x₃ + (_[19])·x₁x₂ + (_[20])·x₁²

julia> @constraint(model, p + q == 1)
(_[1] + _[11] - 1)·1 + (_[2] + _[12])·x₃ + (_[3] + _[13])·x₂ + (_[4] + _[14])·x₁ + (_[5] + _[15])·x₃² + (_[6] + _[16])·x₂x₃ + (_[7] + _[17])·x₂² + (_[8] + _[18])·x₁x₃ + (_[9] + _[19])·x₁x₂ + (_[10] + _[20])·x₁² ∈ PolyJuMP.ZeroPoly()
```

Vectorized constraints can also be used as well as vector of constraints,
named constraints, ...
For instance, if `P` and `Q` are two ``n \times n`` matrices of
polynomials, the following constraints the sum of rows and columns to match:
```julia
@constraint(model, con[i=1:n], sum(P[i, :]) == sum(Q[:, i]))
```
and `con[i]` contains the reference to the constraint between the `i`th row
of `P` and the `i`th column of `Q`.

## Inequality constraints between polynomials

Polynomials can be constrained to be sum-of-squares with the `in` syntax.
For instance, to constrain a polynomial `p` to be sum-of-squares, do
```jldoctest constraint-pq
julia> @constraint(model, p in SOSCone())
(_[1])·1 + (_[2])·x₃ + (_[3])·x₂ + (_[4])·x₁ + (_[5])·x₃² + (_[6])·x₂x₃ + (_[7])·x₂² + (_[8])·x₁x₃ + (_[9])·x₁x₂ + (_[10])·x₁² is SOS
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
```jldoctest constraint-pq
julia> setpolymodule!(model, SumOfSquares)
```
This second approach may be useful if the SumOfSquares JuMP extension need to
be used with another JuMP extension that also has a special model constructor.
A third alternative is the following:
```jldoctest constraint-pq
julia> PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, SOSCone)
SOSCone (alias for NonnegPolyInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle})

julia> PolyJuMP.setdefault!(model, PolyJuMP.PosDefPolyMatrix, SOSMatrixCone)
SOSMatrixCone (alias for PSDMatrixInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle})
```
This approach adds the flexibility to choose the default cone for

* constraints of the form
  `@constraint(mode, ..., some_polynomial ≥ other_polynomial, ...)`
  which is the cone given as default to `PolyJuMP.NonNegPoly`; and
* constraints of the form
  `@constraint(mode, ..., some_matrix_of_polynomial in PSDCone(), ...)`
  or
  `@constraint(mode, ..., some_matrix_of_polynomial >= other_matrix_of_polynomial, PSDCone(), ...)`
  which is the cone given as default to `PolyJuMP.NonNegPolyMatrix`.

For instance, to use the diagonally-dominant-sum-of-squares cone (see
[Ahmadi2017; Definition 2](@cite)) for the first type of contraints, do
```jldoctest constraint-pq
julia> PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, DSOSCone)
DSOSCone (alias for NonnegPolyInnerCone{SumOfSquares.DiagonallyDominantConeTriangle})
```

## Changing the polynomial basis

As introduced in [Choosing a polynomial basis](@ref), there may be numerical
advantages to use another basis than the standard monomial basis when creating
polynomial variables. Similarly, other polynomial bases can be used for
polynomial constraints. However, for constraints, the polynomial space is
determined by the polynomial constrained to be nonnegative. For instance,
consider the constraint:
```jldoctest constraint-xy
julia> using DynamicPolynomials

julia> @polyvar x y
(x, y)

julia> using SumOfSquares

julia> model = SOSModel();

julia> @variable(model, α)
α

julia> @variable(model, β)
β

julia> @constraint(model, α * x^2 + β * y^2 ≥ (α - β) * x * y)
(β)·y² + (-α + β)·xy + (α)·x² is SOS
```
where `α` and `β` are JuMP decision variables and `x` and `y` are polynomial
variables. Since the polynomial is a quadratic form, the sum-of-squares
certificate is also a quadratic form (see [Blekherman2012; Section~3.3.4](@cite)). Hence the
default polynomial basis used for the [Nonnegative polynomial variables]
certificate is `MultivariateBases.SubBasis{MultivariateBases.Monomial}([x, y])`, that is, we search for a positive
semidefinite matrix `Q` such that
```math
\alpha x^2 + \beta y^2 - (\alpha - \beta) x y = X^\top Q X
```
where ``X = (x, y)``.

As the polynomial space is determined by the polynomial being constrained,
only the basis *type* needs to be given. For instance, to use the scaled monomial
basis in the example above, use
```jldoctest constraint-xy
julia> @constraint(model, α * x^2 + β * y^2 ≥ (α - β) * x * y, basis = ScaledMonomial)
(β)·y² + (-α + β)·xy + (α)·x² is SOS
```
In the output, we still see the polynomial in monomial basis.
The transformation to `basis` is done internally, after Newton polytope step.
Although the Newton polytope computation in `Monomial` or `ScaledMonomial` basis
are equivalent,
this is not the case for `Chebyshev` basis for instance.

## Polynomial nonnegativity on a subset of the space

By default, the constraint
```jldoctest constraint-xy
julia> @constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α)
(-α)·1 + (-1)·y² + (2)·xy + (-1)·x² + (1)·y³ + (1)·x³ is SOS
```
constrains the polynomial to be nonnegative for every real numbers `x` and `y`.
However, the set of points `(x, y)` for which the polynomial is constrained
to be nonnegative can be specified by the `domain` keyword:
```jldoctest constraint-xy
julia> S = @set x >= 0 && y >= 0 && x + y >= 1;

julia> @constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α, domain = S)
(-α)·1 + (-1)·y² + (2)·xy + (-1)·x² + (1)·y³ + (1)·x³ is SOS
```
See [this notebook](https://github.com/jump-dev/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb)
for a detailed example.

## Dual of polynomial constraints

The dual of a polynomial constraint `cref` is a moment serie `μ` as defined in
[MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl).
The dual can be obtained with the `dual` function as with classical
dual values in JuMP.
```julia
μ = dual(cref)
```
By dual of a Sum-of-Squares constraint, we may mean different things
and the meaning chosen for `dual` function was chosen for consistency
with the definition of the JuMP `dual` function to ensure that generic
code will work as expected with Sum-of-Squares constraints.
In a Sum-of-Squares constraint, a polynomial $p$ is constraint to
be SOS in some domain defined by polynomial `q_i`.
So `p(x)` is constrained to be equal to
`s(x) = s_0(x) + s_1(x) * q_1(x) + s_2(x) * q_2(x) + ...`
where the `s_i(x)` polynomials are Sum-of-Squares.
The dual of the equality constraint between `p(x)` and `s(x)` is given
by [`SumOfSquares.MultivariateMoments.moments`](@ref).
```julia
μ = moments(cref)
```
Note that the `dual` and `moments` may give different results. For instance,
the output of `dual` only contains the moments corresponding to monomials of `p`
while the output of `moments` may give the  moments of other monomials if `s(x)`
has more monomials than `p(x)`. Besides, if the domain contains polynomial,
equalities, only the  remainder of `p(x) - s(x)` modulo the ideal is constrained
to be zero, see Corollary 2 of [Cox2015](@cite). In that case, the output `moments` is
the dual of the constraint on the remainder so some monomials may have different
moments with `dual` or `moments`.

The dual of the Sum-of-Squares constraint on `s_0(x)`, commonly referred
to as the the matrix of moments can be obtained using [`moment_matrix`](@ref):
```julia
ν = moment_matrix(cref)
```
The `atomic_measure` function of [MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl)
can be used to check if there exists an *atomic* measure (i.e. a measure that is
a sum of Dirac measures) that has the moments given in the the moment matrix
`ν`. This can be used for instance in polynomial optimization (see
[this notebook](https://github.com/jump-dev/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb))
or stability analysis (see
[this notebook](https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/LPJ17e43.ipynb)).

### SAGE extension

To use the SAGE cone in place of the Sum-of-Squares cone for an inequality constraints
between polynomials, use the following:
```julia
import PolyJuMP
PolyJuMP.setpolymodule!(model, PolyJuMP.SAGE)
```
