# Constraints

## Equality constraints between polynomials

Equality between polynomials in
[PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) uses the same syntax as
equality between affine or quadratic expression in
[JuMP](https://github.com/JuliaOpt/JuMP.jl).
For instance, creating two quadratic `n`-variate polynomials `p` and `q` that
must sum up to one can be done as follows:
```jldoctest constraint-pq
julia> n = 3
3

julia> using DynamicPolynomials

julia> @polyvar x[1:n]
(DynamicPolynomials.PolyVar{true}[x₁, x₂, x₃],)

julia> X = monomials(x, 0:2)
10-element MonomialVector{true}:
 x₁²
 x₁x₂
 x₁x₃
 x₂²
 x₂x₃
 x₃²
 x₁
 x₂
 x₃
 1

julia> using SumOfSquares

julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, p, Poly(X))
(noname)x₁² + (noname)x₁x₂ + (noname)x₁x₃ + (noname)x₂² + (noname)x₂x₃ + (noname)x₃² + (noname)x₁ + (noname)x₂ + (noname)x₃ + (noname)

julia> @variable(model, q, Poly(X))
(noname)x₁² + (noname)x₁x₂ + (noname)x₁x₃ + (noname)x₂² + (noname)x₂x₃ + (noname)x₃² + (noname)x₁ + (noname)x₂ + (noname)x₃ + (noname)

julia> @constraint(model, p + q == 1)
(noname + noname)x₁² + (noname + noname)x₁x₂ + (noname + noname)x₁x₃ + (noname + noname)x₂² + (noname + noname)x₂x₃ + (noname + noname)x₃² + (noname + noname)x₁ + (noname + noname)x₂ + (noname + noname)x₃ + (noname + noname - 1) ∈ PolyJuMP.ZeroPoly()
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
(noname)x₁² + (noname)x₁x₂ + (noname)x₁x₃ + (noname)x₂² + (noname)x₂x₃ + (noname)x₃² + (noname)x₁ + (noname)x₂ + (noname)x₃ + (noname) is SOS
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
NonnegPolyInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle}

julia> PolyJuMP.setdefault!(model, PolyJuMP.PosDefPolyMatrix, SOSMatrixCone)
PSDMatrixInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle}
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
```jldoctest constraint-pq
julia> PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, DSOSCone)
NonnegPolyInnerCone{SumOfSquares.DiagonallyDominantConeTriangle}
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

julia> model = SOSModel()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, α)
α

julia> @variable(model, β)
β

julia> @constraint(model, α * x^2 + β * y^2 ≥ (α - β) * x * y)
(α)x² + (-α + β)xy + (β)y² is SOS
```
where `α` and `β` are JuMP decision variables and `x` and `y` are polynomial
variables. Since the polynomial is a quadratic form, the sum-of-squares
certificate is also a quadratic form (see [Section~3.3.4, BPT12]). Hence the
default polynomial basis used for the [Nonnegative polynomial variables]
certificate is `MonomialBasis([x, y])`, that is, we search for a positive
semidefinite matrix `Q` such that
```math
\alpha x^2 + \beta y^2 - (\alpha - \beta) x y = X^\top Q X
```
where ``X = (x, y)``.

As the polynomial space is determined by the polynomial being constrained,
only the basis *type* needs to be given. For instance, to use the scaled monomial
basis, use
```julia
@constraint(model, p ≥ q, basis = ScaledMonomialBasis)
```

## Polynomial nonnegativity on a subset of the space

By default, the constraint
```jldoctest constraint-xy
julia> @constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α)
(1)x³ + (1)y³ + (-1)x² + (2)xy + (-1)y² + (-α) is SOS
```
constrains the polynomial to be nonnegative for every real numbers `x` and `y`.
However, the set of points `(x, y)` for which the polynomial is constrained
to be nonnegative can be specified by the `domain` keyword:
```jldoctest constraint-xy
julia> S = @set x >= 0 && y >= 0 && x + y >= 1;

julia> @constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α, domain = S)
(1)x³ + (1)y³ + (-1)x² + (2)xy + (-1)y² + (-α) is SOS
```
See [this notebook](https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb)
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
by [`SumOfSquares.PolyJuMP.moments`](@ref).
```julia
μ = moments(cref)
```
Note that the `dual` and `moments` may give different results. For instance,
the output of `dual` only contains the moments corresponding to monomials of `p`
while the output of `moments` may give the  moments of other monomials if `s(x)`
has more monomials than `p(x)`. Besides, if the domain contains polynomial,
equalities, only the  remainder of `p(x) - s(x)` modulo the ideal is constrained
to be zero, see Corollary 2 of [CLO13]. In that case, the output `moments` is
the dual of the constraint on the remainder so some monomials may have different
moments with `dual` or `moments`.

The dual of the Sum-of-Squares constraint on `s_0(x)`, commonly referred
to as the the matrix of moments can be obtained using [`moment_matrix`](@ref):
```julia
ν = moment_matrix(cref)
```
The `extractatoms` function of [MultivariateMoments](https://github.com/JuliaAlgebra/MultivariateMoments.jl)
can be used to check if there exists an *atomic* measure (i.e. a measure that is
a sum of Dirac measures) that has the moments given in the the moment matrix
`ν`. This can be used for instance in polynomial optimization (see
[this notebook](https://github.com/JuliaOpt/SumOfSquares.jl/blob/master/examples/Polynomial_Optimization.ipynb))
or stability analysis (see
[this notebook](https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/LPJ17e43.ipynb)).

### References

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.

[CLO13] Cox, D., Little, J., & OShea, D.
*Ideals, varieties, and algorithms: an introduction to computational algebraic geometry and commutative algebra*.
Springer Science & Business Media, **2013**.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*.
ArXiv e-prints, **2017**.

## Reference

Inner approximations of the PSD cone that do not require semidefinite
programming:
```@docs
SumOfSquares.DiagonallyDominantConeTriangle
SumOfSquares.ScaledDiagonallyDominantConeTriangle
```

Approximations of the cone of nonnegative polynomials:
```@docs
SumOfSquares.NonnegPolyInnerCone
SumOfSquares.SOSCone
SumOfSquares.SDSOSCone
SumOfSquares.DSOSCone
```

Approximations of the cone of positive semidefinite polynomial matrices:
```@docs
SumOfSquares.PSDMatrixInnerCone
SumOfSquares.SOSMatrixCone
```

Approximations of the cone of convex polynomials:
```@docs
SumOfSquares.ConvexPolyInnerCone
SumOfSquares.SOSConvexCone
```

Approximations of the cone of copositive matrices:
```@docs
SumOfSquares.CopositiveInner
```

Attributes
```@docs
SumOfSquares.PolyJuMP.MomentsAttribute
SumOfSquares.MultivariateMoments.moments(::SumOfSquares.JuMP.ConstraintRef)
GramMatrix
SumOfSquares.GramMatrixAttribute
gram_matrix
gram_operate
SumOfSquares.MomentMatrixAttribute
moment_matrix
SumOfSquares.CertificateMonomials
certificate_monomials
SumOfSquares.LagrangianMultipliers
lagrangian_multipliers
```

Polynomial basis:
```@docs
SumOfSquares.PolyJuMP.AbstractPolynomialBasis
SumOfSquares.PolyJuMP.MonomialBasis
SumOfSquares.PolyJuMP.ScaledMonomialBasis
SumOfSquares.PolyJuMP.FixedPolynomialBasis
```

Bridges are automatically added using the following utilities:
```@docs
SumOfSquares.PolyJuMP.bridgeable
SumOfSquares.PolyJuMP.bridges
```
