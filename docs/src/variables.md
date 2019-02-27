# Polynomial and Sum-of-Squares variables

## Polynomial variables

While [JuMP](https://github.com/JuliaOpt/JuMP.jl) allows to create decision
variables representing a number whose value needs to be optimized upon by the
optimizer, [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) allows to create
*polynomial* decision variables. In order to do that, we first need to create
polynomial variables with the `@polyvar` macro:
```jldoctest variables
julia> using DynamicPolynomials # or TypedPolynomials, pick your favorite

julia> @polyvar x y
(x, y)
```
Note that these should not be confused with JuMP's decision variables which are
created using the `@variable` macro. The polynomial decision variable that will
be created need to be parametrized by a polynomial basis of finite size.
For instance, if we want to find a quadratic polynomial, we can parametrize it
with all monomials of degree between 0 and 2. Generating a vector with such
monomials can be achieved through the `monomials` function:
```jldoctest variables
julia> X = monomials([x, y], 0:2)
6-element MonomialVector{true}:
 x²
 xy
 y²
 x
 y
 1
```
We can now create our polynomial variable `p` as follows:
```jldoctest variables
julia> using SumOfSquares

julia> model = Model()
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: NO_OPTIMIZER
Solver name: No optimizer attached.

julia> @variable(model, p, Poly(X))
(noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)
```
This creates a vector of decision variables `a` and sets `p` as the scalar
product between `a` and `X`.

Just like with classical JuMP's decision variables, containers of polynomial
variables can be created as follows:
```jldoctest variables
julia> @variable(model, [1:3, 1:4], Poly(X))       # Creates a Matrix
3×4 Array{Polynomial{true,VariableRef},2}:
 (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)  …  (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)
 (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)     (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)
 (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)     (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)

julia> @variable(model, [[:a, :b], -2:2], Poly(X)) # Creates a JuMPArray
2-dimensional DenseAxisArray{DynamicPolynomials.Polynomial{true,VariableRef},2,...} with index sets:
    Dimension 1, Symbol[:a, :b]
    Dimension 2, -2:2
And data, a 2×5 Array{DynamicPolynomials.Polynomial{true,VariableRef},2}:
 (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)  …  (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)
 (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)     (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)

julia> @variable(model, [i=1:3, j=i:3], Poly(X))   # Creates a Dict
JuMP.Containers.SparseAxisArray{Polynomial{true,VariableRef},2,Tuple{Any,Any}} with 6 entries:
  [1, 2]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)
  [2, 3]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)
  [3, 3]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)
  [2, 2]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)
  [1, 1]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)
  [1, 3]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)
```

For more flexibility, polynomials parametrized by decision variables can also
be created "by hand" for instance as follows:
```jldoctest variables
julia> @variable(model, α)
α

julia> @variable(model, β)
β

julia> p = α*x^2 + (α+β)*y^2*x + β*y^3
(α + β)xy² + (β)y³ + (α)x²
```
The coefficients of the polynomial can even be quadratic, e.g.
```jldoctest variables
julia> p = (3α^2+β)*x^2 + (α*β+2β)*y^2*x + β*y^3
(α*β + 2 β)xy² + (β)y³ + (3 α² + β)x²
```
Let me stress again the distinction between `α` and `β` which are *decision*
variables and `x` and `y` which are *polynomial* variables.

## Nonnegative polynomial variables

In order to create a sum-of-squares polynomial variable, the syntax is exactly
the same except `SOSPoly` should be used instead of `Poly`.
For instance, the following code creates a ``3 \times 4`` matrix of
sum-of-squares polynomial variables:
```jldoctest variables
julia> @variable(model, [1:2], SOSPoly(X))
2-element Array{GramMatrix{VariableRef,Monomial{true},MonomialVector{true}},1}:
 GramMatrix{VariableRef,Monomial{true},MonomialVector{true}}(VariableRef[noname noname … noname noname; noname noname … noname noname; … ; noname noname … noname noname; noname noname … noname noname], DynamicPolynomials.Monomial{true}[x², xy, y², x, y, 1])
 GramMatrix{VariableRef,Monomial{true},MonomialVector{true}}(VariableRef[noname noname … noname noname; noname noname … noname noname; … ; noname noname … noname noname; noname noname … noname noname], DynamicPolynomials.Monomial{true}[x², xy, y², x, y, 1])
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

## Choosing a polynomial basis

In the previous section, we show how to create polynomial variables from a
monomial basis. However, the monomial basis is only a particular case of
polynomial basis and while using a different basis of the same space of
polynomial is would give an equivalent program, it might be more stable
numerically (see [Section 3.1.5, BPT12]).

For instance, creating an univariate cubic polynomial variable `p` using the
Chebyshev basis can be done as follows:
```jldoctest variables
julia> cheby_basis = FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])
FixedPolynomialBasis{Polynomial{true,Int64},Array{Polynomial{true,Int64},1}}(DynamicPolynomials.Polynomial{true,Int64}[1, x, 2x² - 1, 4x³ - 3x])

julia> @variable(model, variable_type=Poly(cheby_basis))
(4 noname)x³ + (2 noname)x² + (noname - 3 noname)x + (noname - noname)
```
and to create a quadratic form variable `q` using the *scaled monomial* basis
(see [Section 3.1.5, BPT12]), use the following:
```jldoctest variables
julia> X = monomials([x, y], 2)
3-element MonomialVector{true}:
 x²
 xy
 y²

julia> scaled_basis = ScaledMonomialBasis(X)
ScaledMonomialBasis{Monomial{true},MonomialVector{true}}(DynamicPolynomials.Monomial{true}[x², xy, y²])

julia> @variable(model, variable_type=Poly(scaled_basis))
(noname)x² + (1.4142135623730951 noname)xy + (noname)y²
```
which is equivalent to
```jldoctest variables
julia> scaled_basis = FixedPolynomialBasis([x^2, √2*x*y, y^2])
FixedPolynomialBasis{Term{true,Float64},Array{Term{true,Float64},1}}(DynamicPolynomials.Term{true,Float64}[x², 1.41421xy, y²])

julia> @variable(model, variable_type=Poly(scaled_basis))
(noname)x² + (1.4142135623730951 noname)xy + (noname)y²
```

### References

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.
