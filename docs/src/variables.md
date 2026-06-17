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

julia> model = Model();

julia> @variable(model, p, Poly(X))
(_[1])·1 + (_[2])·y + (_[3])·x + (_[4])·y² + (_[5])·xy + (_[6])·x²
```
This creates a vector of decision variables `a` and sets `p` as the scalar
product between `a` and `X`.

Just like with classical JuMP's decision variables, containers of polynomial
variables can be created as follows:
```jldoctest variables; filter = [r"(Matrix|Vector|DenseAxisArray|SparseAxisArray)\{.*\}" => s"\1{…}"]
julia> @variable(model, [1:3, 1:4], Poly(X))       # Creates a Matrix
3×4 Matrix{…}:
 (_[7])·1 + (_[8])·y + (_[9])·x + (_[10])·y² + (_[11])·xy + (_[12])·x²     …  (_[61])·1 + (_[62])·y + (_[63])·x + (_[64])·y² + (_[65])·xy + (_[66])·x²
 (_[13])·1 + (_[14])·y + (_[15])·x + (_[16])·y² + (_[17])·xy + (_[18])·x²     (_[67])·1 + (_[68])·y + (_[69])·x + (_[70])·y² + (_[71])·xy + (_[72])·x²
 (_[19])·1 + (_[20])·y + (_[21])·x + (_[22])·y² + (_[23])·xy + (_[24])·x²     (_[73])·1 + (_[74])·y + (_[75])·x + (_[76])·y² + (_[77])·xy + (_[78])·x²

julia> @variable(model, [[:a, :b], -2:2], Poly(X)) # Creates a DenseAxisArray
2-dimensional DenseAxisArray{…} with index sets:
    Dimension 1, [:a, :b]
    Dimension 2, -2:2
And data, a 2×5 Matrix{…}:
 (_[79])·1 + (_[80])·y + (_[81])·x + (_[82])·y² + (_[83])·xy + (_[84])·x²  …  (_[127])·1 + (_[128])·y + (_[129])·x + (_[130])·y² + (_[131])·xy + (_[132])·x²
 (_[85])·1 + (_[86])·y + (_[87])·x + (_[88])·y² + (_[89])·xy + (_[90])·x²     (_[133])·1 + (_[134])·y + (_[135])·x + (_[136])·y² + (_[137])·xy + (_[138])·x²

julia> @variable(model, [i=1:3, j=i:3], Poly(X))   # Creates a SparseAxisArray
JuMP.Containers.SparseAxisArray{…} with 6 entries:
  [1, 1]  =  (_[139])·1 + (_[140])·y + (_[141])·x + (_[142])·y^2 + (_[143])·x*y + (_[144])·x^2
  [1, 2]  =  (_[145])·1 + (_[146])·y + (_[147])·x + (_[148])·y^2 + (_[149])·x*y + (_[150])·x^2
  [1, 3]  =  (_[151])·1 + (_[152])·y + (_[153])·x + (_[154])·y^2 + (_[155])·x*y + (_[156])·x^2
  [2, 2]  =  (_[157])·1 + (_[158])·y + (_[159])·x + (_[160])·y^2 + (_[161])·x*y + (_[162])·x^2
  [2, 3]  =  (_[163])·1 + (_[164])·y + (_[165])·x + (_[166])·y^2 + (_[167])·x*y + (_[168])·x^2
  [3, 3]  =  (_[169])·1 + (_[170])·y + (_[171])·x + (_[172])·y^2 + (_[173])·x*y + (_[174])·x^2
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
For instance, the following code creates a vector of 2
sum-of-squares polynomial variables:
```jldoctest variables; filter = [r"(Matrix|Vector|DenseAxisArray|SparseAxisArray)\{.*\}" => s"\1{…}"]
julia> @variable(model, [1:2], SOSPoly(X))
2-element Vector{…}:
 GramMatrix with row/column basis:
 SubBasis{Monomial}([1, y, x, y^2, x*y, x^2])
And entries in a 6×6 SymMatrix{VariableRef}:
 _[177]  _[178]  _[180]  _[183]  _[187]  _[192]
 _[178]  _[179]  _[181]  _[184]  _[188]  _[193]
 _[180]  _[181]  _[182]  _[185]  _[189]  _[194]
 _[183]  _[184]  _[185]  _[186]  _[190]  _[195]
 _[187]  _[188]  _[189]  _[190]  _[191]  _[196]
 _[192]  _[193]  _[194]  _[195]  _[196]  _[197]
 GramMatrix with row/column basis:
 SubBasis{Monomial}([1, y, x, y^2, x*y, x^2])
And entries in a 6×6 SymMatrix{VariableRef}:
 _[198]  _[199]  _[201]  _[204]  _[208]  _[213]
 _[199]  _[200]  _[202]  _[205]  _[209]  _[214]
 _[201]  _[202]  _[203]  _[206]  _[210]  _[215]
 _[204]  _[205]  _[206]  _[207]  _[211]  _[216]
 _[208]  _[209]  _[210]  _[211]  _[212]  _[217]
 _[213]  _[214]  _[215]  _[216]  _[217]  _[218]
```
There is however an *important* difference between the meaning of the
vector of monomials `X` between `Poly` and `SOSPoly`.
In the case of `SOSPoly`, for each SOS polynomial decision variable created,
let's call it `p`, think of it as internally creating a positive semidefinite matrix
`Q` and replacing this `p` by the value of `X' * Q * X` everywhere `p` is used.
So, for instance, if `X` contains all the monomials
of degree 2, then all the monomials of `p` will have degree 4 (i.e. `p` will be
a quartic form).

Similarly, to create diagonally-dominant-sum-of-squares polynomial variables
(see [Ahmadi2017; Definition 3.1](@cite)), use `DSOSPoly(X)`. This creates a diagonally
dominant matrix of variables `Q` and sets the polynomial variables as the value
of `X' * Q * X`.

Finally, to create scaled-diagonally-dominant-sum-of-squares polynomial
variables (see [Ahmadi2017; Definition 3.2](@cite)), use `SDSOSPoly(X)`. This creates a
scaled diagonally dominant matrix of variables `Q` and sets the polynomial
variables as the value of `X' * Q * X`.

## Choosing a polynomial basis

In the previous section, we show how to create polynomial variables from a
monomial basis. However, the monomial basis is only a particular case of
polynomial basis and while using a different basis of the same space of
polynomial would give an equivalent program, it might be more stable
numerically (see [Blekherman2012; Section 3.1.5](@cite)).

For instance, creating an univariate cubic polynomial variable `p` using the
Chebyshev basis can be done as follows:
```jldoctest variables
julia> cheby_basis = SubBasis{Chebyshev}(monomials(x, 0:3))
SubBasis{ChebyshevFirstKind}([1, x, x², x³])

julia> @variable(model, variable_type=Poly(cheby_basis))
(_[219])·ChebyshevFirstKind(1) + (_[220])·ChebyshevFirstKind(x) + (_[221])·ChebyshevFirstKind(x²) + (_[222])·ChebyshevFirstKind(x³)
```
and to create a quadratic form variable `q` using the *scaled monomial* basis
(see [Blekherman2012; Section 3.1.5](@cite)), use the following:
```jldoctest variables
julia> X = monomials([x, y], 2)
3-element MonomialVector{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}:
 y²
 xy
 x²

julia> scaled_basis = SubBasis{ScaledMonomial}(X)
SubBasis{ScaledMonomial}([y², xy, x²])

julia> p = @variable(model, variable_type=Poly(scaled_basis))
(_[223])·ScaledMonomial(y²) + (_[224])·ScaledMonomial(xy) + (_[225])·ScaledMonomial(x²)

julia> polynomial(p)
(_[223])y² + (1.4142135623730951 _[224])xy + (_[225])x²
```
