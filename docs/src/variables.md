# Polynomial and Sum-of-Squares variables

## Polynomial variables

While [JuMP](https://github.com/jump-dev/JuMP.jl) allows to create decision
variables representing a number whose value needs to be optimized upon by the
optimizer, [PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl) allows to create
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
6-element MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}:
 1
 y
 x
 y²
 xy
 x²
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
(_[1]) + (_[2])y + (_[3])x + (_[4])y² + (_[5])xy + (_[6])x²
```
This creates a vector of decision variables `a` and sets `p` as the scalar
product between `a` and `X`.

Just like with classical JuMP's decision variables, containers of polynomial
variables can be created as follows:
```jldoctest variables
julia> @variable(model, [1:3, 1:4], Poly(X))       # Creates a Matrix
3×4 Matrix{Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, VariableRef}}:
 (_[7]) + (_[8])y + (_[9])x + (_[10])y² + (_[11])xy + (_[12])x²     …  (_[61]) + (_[62])y + (_[63])x + (_[64])y² + (_[65])xy + (_[66])x²
 (_[13]) + (_[14])y + (_[15])x + (_[16])y² + (_[17])xy + (_[18])x²     (_[67]) + (_[68])y + (_[69])x + (_[70])y² + (_[71])xy + (_[72])x²
 (_[19]) + (_[20])y + (_[21])x + (_[22])y² + (_[23])xy + (_[24])x²     (_[73]) + (_[74])y + (_[75])x + (_[76])y² + (_[77])xy + (_[78])x²

julia> @variable(model, [[:a, :b], -2:2], Poly(X)) # Creates a DenseAxisArray
2-dimensional DenseAxisArray{DynamicPolynomials.Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, MultivariatePolynomials.Graded{MultivariatePolynomials.LexOrder}, VariableRef},2,...} with index sets:
    Dimension 1, [:a, :b]
    Dimension 2, -2:2
And data, a 2×5 Matrix{DynamicPolynomials.Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, MultivariatePolynomials.Graded{MultivariatePolynomials.LexOrder}, VariableRef}}:
 (_[79]) + (_[80])y + (_[81])x + (_[82])y² + (_[83])xy + (_[84])x²  …  (_[127]) + (_[128])y + (_[129])x + (_[130])y² + (_[131])xy + (_[132])x²
 (_[85]) + (_[86])y + (_[87])x + (_[88])y² + (_[89])xy + (_[90])x²     (_[133]) + (_[134])y + (_[135])x + (_[136])y² + (_[137])xy + (_[138])x²

julia> @variable(model, [i=1:3, j=i:3], Poly(X))   # Creates a SparseAxisArray
JuMP.Containers.SparseAxisArray{Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, VariableRef}, 2, Tuple{Int64, Int64}} with 6 entries:
  [1, 1]  =  (_[139]) + (_[140])*y + (_[141])*x + (_[142])*y^2 + (_[143])*x*y + (_[144])*x^2
  [1, 2]  =  (_[145]) + (_[146])*y + (_[147])*x + (_[148])*y^2 + (_[149])*x*y + (_[150])*x^2
  [1, 3]  =  (_[151]) + (_[152])*y + (_[153])*x + (_[154])*y^2 + (_[155])*x*y + (_[156])*x^2
  [2, 2]  =  (_[157]) + (_[158])*y + (_[159])*x + (_[160])*y^2 + (_[161])*x*y + (_[162])*x^2
  [2, 3]  =  (_[163]) + (_[164])*y + (_[165])*x + (_[166])*y^2 + (_[167])*x*y + (_[168])*x^2
  [3, 3]  =  (_[169]) + (_[170])*y + (_[171])*x + (_[172])*y^2 + (_[173])*x*y + (_[174])*x^2
```

For more flexibility, polynomials parametrized by decision variables can also
be created "by hand" for instance as follows:
```jldoctest variables
julia> @variable(model, α)
α

julia> @variable(model, β)
β

julia> p = α*x^2 + (α+β)*y^2*x + β*y^3
(α)x² + (β)y³ + (α + β)xy²
```
The coefficients of the polynomial can even be quadratic, e.g.
```jldoctest variables
julia> p = (3α^2+β)*x^2 + (α*β+2β)*y^2*x + β*y^3
(3 α² + β)x² + (β)y³ + (α*β + 2 β)xy²
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
2-element Vector{GramMatrix{VariableRef, MonomialBasis{Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}, MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}, AffExpr, SymMatrix{VariableRef}}}:
 GramMatrix{VariableRef, MonomialBasis{Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}, MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}, AffExpr, SymMatrix{VariableRef}}(VariableRef[_[177] _[178] … _[187] _[192]; _[178] _[179] … _[188] _[193]; … ; _[187] _[188] … _[191] _[196]; _[192] _[193] … _[196] _[197]], MonomialBasis{Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}, MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}(Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}[1, y, x, y², xy, x²]))
 GramMatrix{VariableRef, MonomialBasis{Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}, MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}, AffExpr, SymMatrix{VariableRef}}(VariableRef[_[198] _[199] … _[208] _[213]; _[199] _[200] … _[209] _[214]; … ; _[208] _[209] … _[212] _[217]; _[213] _[214] … _[217] _[218]], MonomialBasis{Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}, MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}(Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}[1, y, x, y², xy, x²]))
```
There is however an *important* difference between the meaning of the
vector of monomials `X` between `Poly` and `SOSPoly`. For `SOSPoly`, it
creates a positive semidefinite matrix of variables `Q` and sets `p` as the
value of `X' * Q * X`. That is, for instance, if `X` contains all the monomials
of degree 2, then all the monomials of `p` will have degree 4 (i.e. `p` will be
a quartic form).

Similarly, to create diagonally-dominant-sum-of-squares polynomial variables
(see [Definition 3.1, AM17]), use `DSOSPoly(X)`. This creates a diagonally
dominant matrix of variables `Q` and sets the polynomial variables as the value
of `X' * Q * X`.

Finally, to create scaled-diagonally-dominant-sum-of-squares polynomial
variables (see [Definition 3.2, AM17]), use `SDSOSPoly(X)`. This creates a
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
FixedPolynomialBasis{Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}, Vector{Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}}}(Polynomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}, Int64}[1, x, -1 + 2x², -3x + 4x³])

julia> @variable(model, variable_type=Poly(cheby_basis))
(_[219] - _[221]) + (_[220] - 3 _[222])x + (2 _[221])x² + (4 _[222])x³
```
and to create a quadratic form variable `q` using the *scaled monomial* basis
(see [Section 3.1.5, BPT12]), use the following:
```jldoctest variables
julia> X = monomials([x, y], 2)
3-element MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}:
 y²
 xy
 x²

julia> scaled_basis = ScaledMonomialBasis(X)
ScaledMonomialBasis{Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}, MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}(Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}[y², xy, x²])

julia> @variable(model, variable_type=Poly(scaled_basis))
(_[223])y² + (1.4142135623730951 _[224])xy + (_[225])x²
```
which is equivalent to
```jldoctest variables
julia> scaled_basis = FixedPolynomialBasis([x^2, √2*x*y, y^2])
FixedPolynomialBasis{Term{Float64, Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}, Vector{Term{Float64, Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}}}(Term{Float64, Monomial{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}}[x², 1.4142135623730951xy, y²])

julia> @variable(model, variable_type=Poly(scaled_basis))
(_[228])y² + (1.4142135623730951 _[227])xy + (_[226])x²
```

### References

[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R.
*Semidefinite Optimization and Convex Algebraic Geometry*.
Society for Industrial and Applied Mathematics, **2012**.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.

## Reference

```@docs
SumOfSquares.PolyJuMP.Poly
```

Variable bridges:
```@docs
SumOfSquares.Bridges.Variable.ScaledDiagonallyDominantBridge
```
