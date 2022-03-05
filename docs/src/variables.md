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
(_[1])x² + (_[2])xy + (_[3])y² + (_[4])x + (_[5])y + (_[6])
```
This creates a vector of decision variables `a` and sets `p` as the scalar
product between `a` and `X`.

Just like with classical JuMP's decision variables, containers of polynomial
variables can be created as follows:
```jldoctest variables
julia> @variable(model, [1:3, 1:4], Poly(X))       # Creates a Matrix
3×4 Matrix{Polynomial{true, VariableRef}}:
 (_[7])x² + (_[8])xy + (_[9])y² + (_[10])x + (_[11])y + (_[12])     …  (_[61])x² + (_[62])xy + (_[63])y² + (_[64])x + (_[65])y + (_[66])
 (_[13])x² + (_[14])xy + (_[15])y² + (_[16])x + (_[17])y + (_[18])     (_[67])x² + (_[68])xy + (_[69])y² + (_[70])x + (_[71])y + (_[72])
 (_[19])x² + (_[20])xy + (_[21])y² + (_[22])x + (_[23])y + (_[24])     (_[73])x² + (_[74])xy + (_[75])y² + (_[76])x + (_[77])y + (_[78])

julia> @variable(model, [[:a, :b], -2:2], Poly(X)) # Creates a DenseAxisArray
2-dimensional DenseAxisArray{DynamicPolynomials.Polynomial{true, VariableRef},2,...} with index sets:
    Dimension 1, [:a, :b]
    Dimension 2, -2:2
And data, a 2×5 Matrix{DynamicPolynomials.Polynomial{true, VariableRef}}:
 (_[79])x² + (_[80])xy + (_[81])y² + (_[82])x + (_[83])y + (_[84])  …  (_[127])x² + (_[128])xy + (_[129])y² + (_[130])x + (_[131])y + (_[132])
 (_[85])x² + (_[86])xy + (_[87])y² + (_[88])x + (_[89])y + (_[90])     (_[133])x² + (_[134])xy + (_[135])y² + (_[136])x + (_[137])y + (_[138])

julia> @variable(model, [i=1:3, j=i:3], Poly(X))   # Creates a SparseAxisArray
JuMP.Containers.SparseAxisArray{Polynomial{true, VariableRef}, 2, Tuple{Int64, Int64}} with 6 entries:
  [1, 1]  =  (_[139])*x^2 + (_[140])*x*y + (_[141])*y^2 + (_[142])*x + (_[143])*y + (_[144])
  [1, 2]  =  (_[145])*x^2 + (_[146])*x*y + (_[147])*y^2 + (_[148])*x + (_[149])*y + (_[150])
  [1, 3]  =  (_[151])*x^2 + (_[152])*x*y + (_[153])*y^2 + (_[154])*x + (_[155])*y + (_[156])
  [2, 2]  =  (_[157])*x^2 + (_[158])*x*y + (_[159])*y^2 + (_[160])*x + (_[161])*y + (_[162])
  [2, 3]  =  (_[163])*x^2 + (_[164])*x*y + (_[165])*y^2 + (_[166])*x + (_[167])*y + (_[168])
  [3, 3]  =  (_[169])*x^2 + (_[170])*x*y + (_[171])*y^2 + (_[172])*x + (_[173])*y + (_[174])
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
2-element Vector{GramMatrix{VariableRef, MonomialBasis{Monomial{true}, MonomialVector{true}}, AffExpr, SymMatrix{VariableRef}}}:
 GramMatrix{VariableRef, MonomialBasis{Monomial{true}, MonomialVector{true}}, AffExpr, SymMatrix{VariableRef}}(VariableRef[_[177] _[178] … _[187] _[192]; _[178] _[179] … _[188] _[193]; … ; _[187] _[188] … _[191] _[196]; _[192] _[193] … _[196] _[197]], MonomialBasis{Monomial{true}, MonomialVector{true}}(DynamicPolynomials.Monomial{true}[x², xy, y², x, y, 1]))
 GramMatrix{VariableRef, MonomialBasis{Monomial{true}, MonomialVector{true}}, AffExpr, SymMatrix{VariableRef}}(VariableRef[_[198] _[199] … _[208] _[213]; _[199] _[200] … _[209] _[214]; … ; _[208] _[209] … _[212] _[217]; _[213] _[214] … _[217] _[218]], MonomialBasis{Monomial{true}, MonomialVector{true}}(DynamicPolynomials.Monomial{true}[x², xy, y², x, y, 1]))
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
FixedPolynomialBasis{Polynomial{true, Int64}, Vector{Polynomial{true, Int64}}}(DynamicPolynomials.Polynomial{true, Int64}[1, x, 2x² - 1, 4x³ - 3x])

julia> @variable(model, variable_type=Poly(cheby_basis))
(4 _[222])x³ + (2 _[221])x² + (_[220] - 3 _[222])x + (_[219] - _[221])
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
ScaledMonomialBasis{Monomial{true}, MonomialVector{true}}(DynamicPolynomials.Monomial{true}[x², xy, y²])

julia> @variable(model, variable_type=Poly(scaled_basis))
(_[223])x² + (1.4142135623730951 _[224])xy + (_[225])y²
```
which is equivalent to
```jldoctest variables
julia> scaled_basis = FixedPolynomialBasis([x^2, √2*x*y, y^2])
FixedPolynomialBasis{Term{true, Float64}, Vector{Term{true, Float64}}}(DynamicPolynomials.Term{true, Float64}[x², 1.4142135623730951xy, y²])

julia> @variable(model, variable_type=Poly(scaled_basis))
(_[226])x² + (1.4142135623730951 _[227])xy + (_[228])y²
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
