var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "#SumOfSquares-–-Sum-of-Squares-Programming-for-Julia-1",
    "page": "Index",
    "title": "SumOfSquares –- Sum of Squares Programming for Julia",
    "category": "section",
    "text": "SumOfSquares implements Sum of Squares reformulation for PolyJuMP, enabling the creation of sum of squares variables and constraints in JuMP.The polynomial can be represented by any library implementing the MultivariatePolynomial.jl interface. That is, you can currently choose between DynamicPolynomials and TypedPolynomials. As a rule of thumb, if you know at compile time (or at the time you are writing your code), the number of variable and that this number is small, use TypedPolynomials, otherwise, use DynamicPolynomials.Some presentations on, or using, SumOfSquares:Benoît Legat at the JuMP Meetup 2017 [Slides] [Video]\nJoey Huchette at SIAM Opt 2017The following example shows how to find lower bounds for the Goldstein-Price function using this package.using SumOfSquares\nusing DynamicPolynomials\nusing MosekTools\n\n# Create symbolic variables (not JuMP decision variables)\n@polyvar x1 x2\n\n# Create a Sum of Squares JuMP model with the Mosek solver\nmodel = SOSModel(with_optimizer(Mosek.Optimizer))\n\n# Create a JuMP decision variable for the lower bound\n@variable(model, γ)\n\n# f(x) is the Goldstein-Price function\nf1 = x1+x2+1\nf2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2\nf3 = 2*x1-3*x2\nf4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2\n\nf = (1+f1^2*f2)*(30+f3^2*f4)\n\n# Constraints f(x) - γ to be sum of squares\n@constraint(model, f >= γ)\n\n@objective(model, Max, γ)\n\noptimize!(model)\n\n# The lower bound found is 3\nprintln(objective_value(model))"
},

{
    "location": "#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"sumofsquares.md\", \"variables.md\", \"constraints.md\"]\nDepth = 2"
},

{
    "location": "sumofsquares/#",
    "page": "Sum-of-Squares Programming",
    "title": "Sum-of-Squares Programming",
    "category": "page",
    "text": ""
},

{
    "location": "sumofsquares/#Sum-of-Squares-Programming-1",
    "page": "Sum-of-Squares Programming",
    "title": "Sum-of-Squares Programming",
    "category": "section",
    "text": "This section contains a brief introduction to Sum-of-Squares Programming. For more details, see [BPT12, Las09, Lau09]."
},

{
    "location": "sumofsquares/#Quadratic-forms-and-Semidefinite-programming-1",
    "page": "Sum-of-Squares Programming",
    "title": "Quadratic forms and Semidefinite programming",
    "category": "section",
    "text": "The positive semidefiniteness of a n times n real symmetric matrix Q is equivalent to the nonnegativity of the quadratic form p(x) = x^top Q x for all vector x in mathbbR^n. For instance, the polynomialx_1^2 + 2x_1x_2 + 5x_2^2 + 4x_2x_3 + x_3^2 = x^top beginpmatrix1  1  01  5  2 0  2  1endpmatrix xis nonnegative since the matrix of the right-hand side is positive semidefinite. Moreover, a certificate of nonnegativity can be extracted from the Cholesky decomposition of the matrix:(x_1 + x_2)^2 + (2x_2 + x_3)^2 = x^top beginpmatrix1  1  00  2  1endpmatrix^top beginpmatrix1  1  00  2  1endpmatrix x"
},

{
    "location": "sumofsquares/#Polynomial-nonnegativity-and-Semidefinite-programming-1",
    "page": "Sum-of-Squares Programming",
    "title": "Polynomial nonnegativity and Semidefinite programming",
    "category": "section",
    "text": "This can be generalized to a polynomial of arbitrary degree. A polynomial p(x) is nonnegative is it can be rewritten as p(x) = X^top Q X where Q is a real symmetric positive semidefinite matrix and X is a vector of monomials.For instancex_1^2 + 2x_1^2x_2 + 5x_1^2x_2^2 + 4x_1x_2^2 + x_2^2 = X^top beginpmatrix1  1  01  5  2 0  2  1endpmatrix Xwhere X = (x_1 x_1x_2 x_2) Similarly to the previous section, the Cholesky factorization of the matrix can be used to extract a sum of squares certificate of nonnegativity for the polynomial:(x_1 + x_1x_2)^2 + (2x_1x_2 + x_2)^2 = X^top beginpmatrix1  1  00  2  1endpmatrix^top beginpmatrix1  1  00  2  1endpmatrix X"
},

{
    "location": "sumofsquares/#When-is-nonnegativity-equivalent-to-sum-of-squares-?-1",
    "page": "Sum-of-Squares Programming",
    "title": "When is nonnegativity equivalent to sum of squares ?",
    "category": "section",
    "text": "Determining whether a polynomial is nonnegative is NP-hard. The condition of the previous section was only sufficient, that is, there exists nonnegative polynomials that are not sums of squares. Hilbert showed in 1888 that there are exactly 3 cases for which there is equivalence between the nonnegativity of the polynomials of n variables and degree 2d and the existence of a sum of squares decomposition.n = 1 : Univariate polynomials\n2d = 2 : Quadratic polynomials\nn = 2, 2d = 4 : Bivariate quarticsThe first explicit example of polynomial that was not a sum of squares was given by Motzkin in 1967:x_1^4x_2^2 + x_1^2x_2^4 + 1 - 3x_1^2x_2^2 geq 0 quad forall xWhile it is not a sum of squares, it can still be certified to be nonnegative using sum-of-squares programming by identifying it with a rational sum-of-squares decomposition. These facts can be verified numerically using this package as detailed in the motzkin notebook."
},

{
    "location": "sumofsquares/#References-1",
    "page": "Sum-of-Squares Programming",
    "title": "References",
    "category": "section",
    "text": "[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.[Las09] Lasserre, J. B. Moments, positive polynomials and their applications World Scientific, 2009.[Lau09] Laurent, M. Sums of squares, moment matrices and optimization over polynomials Emerging applications of algebraic geometry, Springer, 2009, 157-270."
},

{
    "location": "variables/#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": ""
},

{
    "location": "variables/#Polynomial-and-Sum-of-Squares-variables-1",
    "page": "Variables",
    "title": "Polynomial and Sum-of-Squares variables",
    "category": "section",
    "text": ""
},

{
    "location": "variables/#Polynomial-variables-1",
    "page": "Variables",
    "title": "Polynomial variables",
    "category": "section",
    "text": "While JuMP allows to create decision variables representing a number whose value needs to be optimized upon by the optimizer, PolyJuMP allows to create polynomial decision variables. In order to do that, we first need to create polynomial variables with the @polyvar macro:julia> using DynamicPolynomials # or TypedPolynomials, pick your favorite\n\njulia> @polyvar x y\n(x, y)Note that these should not be confused with JuMP\'s decision variables which are created using the @variable macro. The polynomial decision variable that will be created need to be parametrized by a polynomial basis of finite size. For instance, if we want to find a quadratic polynomial, we can parametrize it with all monomials of degree between 0 and 2. Generating a vector with such monomials can be achieved through the monomials function:julia> X = monomials([x, y], 0:2)\n6-element MonomialVector{true}:\n x²\n xy\n y²\n x\n y\n 1We can now create our polynomial variable p as follows:julia> using SumOfSquares\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, p, Poly(X))\n(noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)This creates a vector of decision variables a and sets p as the scalar product between a and X.Just like with classical JuMP\'s decision variables, containers of polynomial variables can be created as follows:julia> @variable(model, [1:3, 1:4], Poly(X))       # Creates a Matrix\n3×4 Array{Polynomial{true,VariableRef},2}:\n (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)  …  (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)\n (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)     (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)\n (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)     (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)\n\njulia> @variable(model, [[:a, :b], -2:2], Poly(X)) # Creates a JuMPArray\n2-dimensional DenseAxisArray{DynamicPolynomials.Polynomial{true,VariableRef},2,...} with index sets:\n    Dimension 1, Symbol[:a, :b]\n    Dimension 2, -2:2\nAnd data, a 2×5 Array{DynamicPolynomials.Polynomial{true,VariableRef},2}:\n (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)  …  (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)\n (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)     (noname)x² + (noname)xy + (noname)y² + (noname)x + (noname)y + (noname)\n\njulia> @variable(model, [i=1:3, j=i:3], Poly(X))   # Creates a Dict\nJuMP.Containers.SparseAxisArray{Polynomial{true,VariableRef},2,Tuple{Any,Any}} with 6 entries:\n  [1, 2]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)\n  [2, 3]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)\n  [3, 3]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)\n  [2, 2]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)\n  [1, 1]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)\n  [1, 3]  =  (noname)*x^2 + (noname)*x*y + (noname)*y^2 + (noname)*x + (noname)*y + (noname)For more flexibility, polynomials parametrized by decision variables can also be created \"by hand\" for instance as follows:julia> @variable(model, α)\nα\n\njulia> @variable(model, β)\nβ\n\njulia> p = α*x^2 + (α+β)*y^2*x + β*y^3\n(α + β)xy² + (β)y³ + (α)x²The coefficients of the polynomial can even be quadratic, e.g.julia> p = (3α^2+β)*x^2 + (α*β+2β)*y^2*x + β*y^3\n(α*β + 2 β)xy² + (β)y³ + (3 α² + β)x²Let me stress again the distinction between α and β which are decision variables and x and y which are polynomial variables."
},

{
    "location": "variables/#Nonnegative-polynomial-variables-1",
    "page": "Variables",
    "title": "Nonnegative polynomial variables",
    "category": "section",
    "text": "In order to create a sum-of-squares polynomial variable, the syntax is exactly the same except SOSPoly should be used instead of Poly. For instance, the following code creates a 3 times 4 matrix of sum-of-squares polynomial variables:julia> @variable(model, [1:2], SOSPoly(X))\n2-element Array{GramMatrix{VariableRef,Monomial{true},MonomialVector{true}},1}:\n GramMatrix{VariableRef,Monomial{true},MonomialVector{true}}(VariableRef[noname noname … noname noname; noname noname … noname noname; … ; noname noname … noname noname; noname noname … noname noname], DynamicPolynomials.Monomial{true}[x², xy, y², x, y, 1])\n GramMatrix{VariableRef,Monomial{true},MonomialVector{true}}(VariableRef[noname noname … noname noname; noname noname … noname noname; … ; noname noname … noname noname; noname noname … noname noname], DynamicPolynomials.Monomial{true}[x², xy, y², x, y, 1])There is however an important difference between the signification of the vector of monomials X between Poly and SOSPoly. For SOSPoly, it creates a positive semidefinite matrix of variables Q and sets p as the value of X\' * Q * X. That is, for instance, if X contains all the monomials of degree 2, then all the monomials of p will have degree 4 (i.e. p will be a quartic form).Similarly, to create diagonally-dominant-sum-of-squares polynomial variables (see [Definition 2, AM17]), use DSOSPoly(X). This creates a diagonally dominant matrix of variables Q and sets the polynomial variables as the value of X\' * Q * X.Finally, to create scaled-diagonally-dominant-sum-of-squares polynomial variables (see [Definition 2, AM17]), use DSOSPoly(X). This creates a scaled diagonally dominant matrix of variables Q and sets the polynomial variables as the value of X\' * Q * X."
},

{
    "location": "variables/#Choosing-a-polynomial-basis-1",
    "page": "Variables",
    "title": "Choosing a polynomial basis",
    "category": "section",
    "text": "In the previous section, we show how to create polynomial variables from a monomial basis. However, the monomial basis is only a particular case of polynomial basis and while using a different basis of the same space of polynomial is would give an equivalent program, it might be more stable numerically (see [Section 3.1.5, BPT12]).For instance, creating an univariate cubic polynomial variable p using the Chebyshev basis can be done as follows:julia> cheby_basis = FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])\nFixedPolynomialBasis{Polynomial{true,Int64},Array{Polynomial{true,Int64},1}}(DynamicPolynomials.Polynomial{true,Int64}[1, x, 2x² - 1, 4x³ - 3x])\n\njulia> @variable(model, variable_type=Poly(cheby_basis))\n(4 noname)x³ + (2 noname)x² + (noname - 3 noname)x + (noname - noname)and to create a quadratic form variable q using the scaled monomial basis (see [Section 3.1.5, BPT12]), use the following:julia> X = monomials([x, y], 2)\n3-element MonomialVector{true}:\n x²\n xy\n y²\n\njulia> scaled_basis = ScaledMonomialBasis(X)\nScaledMonomialBasis{Monomial{true},MonomialVector{true}}(DynamicPolynomials.Monomial{true}[x², xy, y²])\n\njulia> @variable(model, variable_type=Poly(scaled_basis))\n(noname)x² + (1.4142135623730951 noname)xy + (noname)y²which is equivalent tojulia> scaled_basis = FixedPolynomialBasis([x^2, √2*x*y, y^2])\nFixedPolynomialBasis{Term{true,Float64},Array{Term{true,Float64},1}}(DynamicPolynomials.Term{true,Float64}[x², 1.41421xy, y²])\n\njulia> @variable(model, variable_type=Poly(scaled_basis))\n(noname)x² + (1.4142135623730951 noname)xy + (noname)y²"
},

{
    "location": "variables/#References-1",
    "page": "Variables",
    "title": "References",
    "category": "section",
    "text": "[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017."
},

{
    "location": "constraints/#",
    "page": "Constraints",
    "title": "Constraints",
    "category": "page",
    "text": ""
},

{
    "location": "constraints/#Constraints-1",
    "page": "Constraints",
    "title": "Constraints",
    "category": "section",
    "text": ""
},

{
    "location": "constraints/#Equality-constraints-between-polynomials-1",
    "page": "Constraints",
    "title": "Equality constraints between polynomials",
    "category": "section",
    "text": "Equality between polynomials in PolyJuMP uses the same syntax as equality between affine or quadratic expression in JuMP. For instance, creating two quadratic n-variate polynomials p and q that must sum up to one can be done as follows:julia> n = 3\n3\n\njulia> using DynamicPolynomials\n\njulia> @polyvar x[1:n]\n(DynamicPolynomials.PolyVar{true}[x₁, x₂, x₃],)\n\njulia> X = monomials(x, 0:2)\n10-element MonomialVector{true}:\n x₁²\n x₁x₂\n x₁x₃\n x₂²\n x₂x₃\n x₃²\n x₁\n x₂\n x₃\n 1\n\njulia> using SumOfSquares\n\njulia> model = Model()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, p, Poly(X))\n(noname)x₁² + (noname)x₁x₂ + (noname)x₁x₃ + (noname)x₂² + (noname)x₂x₃ + (noname)x₃² + (noname)x₁ + (noname)x₂ + (noname)x₃ + (noname)\n\njulia> @variable(model, q, Poly(X))\n(noname)x₁² + (noname)x₁x₂ + (noname)x₁x₃ + (noname)x₂² + (noname)x₂x₃ + (noname)x₃² + (noname)x₁ + (noname)x₂ + (noname)x₃ + (noname)\n\njulia> @constraint(model, p + q == 1)\n(noname + noname)x₁² + (noname + noname)x₁x₂ + (noname + noname)x₁x₃ + (noname + noname)x₂² + (noname + noname)x₂x₃ + (noname + noname)x₃² + (noname + noname)x₁ + (noname + noname)x₂ + (noname + noname)x₃ + (noname + noname - 1) ∈ PolyJuMP.ZeroPoly()Vectorized constraints can also be used as well as vector of constraints, named constraints, ... For instance, if P and Q are two n times n matrices of polynomials, the following constraints the sum of rows and columns to match:@constraint(model, con[i=1:n], sum(P[i, :]) == sum(Q[:, i]))and con[i] contains the reference to the constraint between the ith row of P and the ith column of Q."
},

{
    "location": "constraints/#Inequality-constraints-between-polynomials-1",
    "page": "Constraints",
    "title": "Inequality constraints between polynomials",
    "category": "section",
    "text": "Polynomials can be constrained to be sum-of-squares with the in syntax. For instance, to constrain a polynomial p to be sum-of-squares, dojulia> @constraint(model, p in SOSCone())\n(noname)x₁² + (noname)x₁x₂ + (noname)x₁x₃ + (noname)x₂² + (noname)x₂x₃ + (noname)x₃² + (noname)x₁ + (noname)x₂ + (noname)x₃ + (noname) is SOS"
},

{
    "location": "constraints/#Automatically-interpreting-polynomial-nonnegativity-as-a-sum-of-squares-constraint-1",
    "page": "Constraints",
    "title": "Automatically interpreting polynomial nonnegativity as a sum-of-squares constraint",
    "category": "section",
    "text": "As detailed in When is nonnegativity equivalent to sum of squares ?, the nonnegativity of a polynomial is not equivalent to the existence of a sum-of-squares decomposition. However, if explicitely specified, nonnegativity constraints can be automatically interpreted as sum-of-squares constraints. The simplest way to do that is to create the model withmodel = SOSModel(...)instead ofmodel = Model(...)An alternative equivalent way is to call setpolymodule! after creating the model:julia> setpolymodule!(model, SumOfSquares)This second approach may be useful if the SumOfSquares JuMP extension need to be used with another JuMP extension that also has a special model constructor. A third alternative is the following:julia> PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, SOSCone)\nNonnegPolyInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle}\n\njulia> PolyJuMP.setdefault!(model, PolyJuMP.PosDefPolyMatrix, SOSMatrixCone)\nPSDMatrixInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle}This approach adds the flexibility to choose the default cone forconstraints of the form @constraint(mode, ..., some_polynomial ≥ other_polynomial, ...) which is the cone given as default to PolyJuMP.NonNegPoly; and\nconstraints of the form @constraint(mode, ..., some_matrix_of_polynomial in PSDCone(), ...) or @SDconstraint(mode, ..., some_matrix_of_polynomial ⪰ other_matrix_of_polynomial, ...) which is the cone given as default to PolyJuMP.NonNegPolyMatrix.For instance, to use the diagonally-dominant-sum-of-squares cone (see [Definition 2, AM17]) for the first type of contraints, dojulia> PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, DSOSCone)\nNonnegPolyInnerCone{SumOfSquares.DiagonallyDominantConeTriangle}"
},

{
    "location": "constraints/#Changing-the-polynomial-basis-1",
    "page": "Constraints",
    "title": "Changing the polynomial basis",
    "category": "section",
    "text": "As introduced in Choosing a polynomial basis, there may be numerical advantages to use another basis than the standard monomial basis when creating polynomial variables. Similarly, other polynomial bases can be used for polynomial constraints. However, for constraints, the polynomial space is determined by the polynomial constrained to be nonnegative. For instance, consider the constraint:julia> using DynamicPolynomials\n\njulia> @polyvar x y\n(x, y)\n\njulia> using SumOfSquares\n\njulia> model = SOSModel()\nA JuMP Model\nFeasibility problem with:\nVariables: 0\nModel mode: AUTOMATIC\nCachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer attached.\n\njulia> @variable(model, α)\nα\n\njulia> @variable(model, β)\nβ\n\njulia> @constraint(model, α * x^2 + β * y^2 ≥ (α - β) * x * y)\n(α)x² + (-α + β)xy + (β)y² is SOSwhere α and β are JuMP decision variables and x and y are polynomial variables. Since the polynomial is a quadratic form, the sum-of-squares certificate is also a quadratic form (see [Section~3.3.4, BPT12]). Hence the default polynomial basis used for the [Nonnegative polynomial variables] certificate is MonomialBasis([x, y]), that is, we search for a positive semidefinite matrix Q such thatalpha x^2 + beta y^2 - (alpha - beta) x y = X^top Q Xwhere X = (x y).As the polynomial space is determined by the polynomial being constrained, only the basis type needs to be given. For instance, to use the scaled monomial basis, use@constraint(model, p ≥ q, basis = ScaledMonomialBasis)"
},

{
    "location": "constraints/#Polynomial-nonnegativity-on-a-subset-of-the-space-1",
    "page": "Constraints",
    "title": "Polynomial nonnegativity on a subset of the space",
    "category": "section",
    "text": "By default, the constraintjulia> @constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α)\n(1)x³ + (1)y³ + (-1)x² + (2)xy + (-1)y² + (-α) is SOSconstrains the polynomial to be nonnegative for every real numbers x and y. However, the set of points (x, y) for which the polynomial is constrained to be nonnegative can be specified by the domain keyword:julia> S = @set x >= 0 && y >= 0 && x + y >= 1;\n\njulia> @constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α, domain = S)\n(1)x³ + (1)y³ + (-1)x² + (2)xy + (-1)y² + (-α) is SOSSee this notebook for a detailed example."
},

{
    "location": "constraints/#Dual-of-polynomial-constraints-1",
    "page": "Constraints",
    "title": "Dual of polynomial constraints",
    "category": "section",
    "text": "The dual of a polynomial constraint cref is a moment serie μ as defined in MultivariateMoments. The dual can be obtained with the dual function as with classical dual values in JuMP. The matrix of moments can be obtained using moment_matrix:μ = dual(cref)\nν = moment_matrix(cref)The extractatoms function of MultivariateMoments can be used to check if there exists an atomic measure (i.e. a measure that is a sum of Dirac measures) that has the moments given in ν. This can be used for instance in polynomial optimization (see this notebook) or stability analysis (see this notebook)."
},

{
    "location": "constraints/#References-1",
    "page": "Constraints",
    "title": "References",
    "category": "section",
    "text": "[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017."
},

{
    "location": "constraints/#SumOfSquares.DiagonallyDominantConeTriangle",
    "page": "Constraints",
    "title": "SumOfSquares.DiagonallyDominantConeTriangle",
    "category": "type",
    "text": "struct DiagonallyDominantConeTriangle <: MatrixConeTriangle\n    side_dimension::Int\nend\n\nSee Definition 4 of [AM17] for a precise definition of the last two items.\n\n[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.ScaledDiagonallyDominantConeTriangle",
    "page": "Constraints",
    "title": "SumOfSquares.ScaledDiagonallyDominantConeTriangle",
    "category": "type",
    "text": "struct ScaledDiagonallyDominantConeTriangle <: MatrixConeTriangle\n    side_dimension::Int\nend\n\nSee Definition 4 of [AM17] for a precise definition of the last two items.\n\n[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.NonnegPolyInnerCone",
    "page": "Constraints",
    "title": "SumOfSquares.NonnegPolyInnerCone",
    "category": "type",
    "text": "struct NonnegPolyInnerCone{MCT <: MOI.AbstractVectorSet}\nend\n\nInner approximation of the cone of nonnegative polynomials defined by the set of polynomials of the form\n\nX^\\top Q X\n\nwhere X is a vector of monomials and Q is a symmetric matrix that belongs to the cone psd_inner.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.SOSCone",
    "page": "Constraints",
    "title": "SumOfSquares.SOSCone",
    "category": "type",
    "text": "const SOSCone = NonnegPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}\n\nSum-of-squares cone; see NonnegPolyInnerCone.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.SDSOSCone",
    "page": "Constraints",
    "title": "SumOfSquares.SDSOSCone",
    "category": "type",
    "text": "const SDSOSCone = NonnegPolyInnerCone{ScaledDiagonallyDominantConeTriangle}\n\nScaled-diagonally-dominant-sum-of-squares cone; see [Definition 2, AM17] and NonnegPolyInnerCone.\n\n[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.DSOSCone",
    "page": "Constraints",
    "title": "SumOfSquares.DSOSCone",
    "category": "type",
    "text": "const DSOSCone = NonnegPolyInnerCone{DiagonallyDominantConeTriangle}\n\nDiagonally-dominant-sum-of-squares cone; see [Definition 2, AM17] and NonnegPolyInnerCone.\n\n[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.PSDMatrixInnerCone",
    "page": "Constraints",
    "title": "SumOfSquares.PSDMatrixInnerCone",
    "category": "type",
    "text": "struct PSDMatrixInnerCone{MCT <: MOI.AbstractVectorSet} <: PolyJuMP.PolynomialSet\nend\n\nInner approximation of the cone of polynomial matrices P(x) that are positive semidefinite for any x defined by the set of n times n polynomial matrices such that the polynomial y^top P(x) y belongs to NonnegPolyInnerCone{MCT} where y is a vector of n auxiliary polynomial variables.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.SOSMatrixCone",
    "page": "Constraints",
    "title": "SumOfSquares.SOSMatrixCone",
    "category": "type",
    "text": "const SOSMatrixCone = PSDMatrixInnerCone{MOI.PositiveSemidefiniteConeTriangle}\n\nSum-of-squares matrices cone; see [Section 3.3.2, BPT12] and PSDMatrixInnerCone.\n\n[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.ConvexPolyInnerCone",
    "page": "Constraints",
    "title": "SumOfSquares.ConvexPolyInnerCone",
    "category": "type",
    "text": "struct ConvexPolyInnerCone{MCT} <: PolyJuMP.PolynomialSet end\n\nInner approximation of the set of convex polynomials defined by the set of polynomials such that their hessian belongs to to the set PSDMatrixInnerCone{MCT}()\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.SOSConvexCone",
    "page": "Constraints",
    "title": "SumOfSquares.SOSConvexCone",
    "category": "type",
    "text": "const SOSConvexCone = ConvexPolyInnerCone{MOI.PositiveSemidefiniteConeTriangle}\n\nSum-of-squares convex polynomials cone; see [Section 3.3.3, BPT12] and ConvexPolyInnerCone.\n\n[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.CopositiveInner",
    "page": "Constraints",
    "title": "SumOfSquares.CopositiveInner",
    "category": "type",
    "text": "struct CopositiveInner{S} <: PolyJuMP.PolynomialSet\n    # Inner approximation of the PSD cone, i.e. typically either\n    # `SOSCone`, `DSOSCone` or `SDSOSCone`,\n    psd_inner::S\nend\n\nA symmetric matrix Q is copositive if x^top Q x ge 0 for all vector x in the nonnegative orthant. Checking copositivity is a co-NP-complete problem [MK87] and this cone is only the inner approximation of the cone of copositive symmetric matrices given by Minknowski sum of psd_inner and the cone of symmetric matrices with nonnegative entries (the diagonal entries can be chosen to be zero) [Lemma 3.164, BPT12].\n\nThe matrix with nonnegative entries can be interpreted as lagrangian multipliers. For instance,\n\n@polyvar x y\n@constraint(model, x^2 - 2x*y + y^2 in CopositiveInner(SOSCone()))\n\nis equivalent to\n\n# Matrix that we require to be copositive\nQ = [ 1 -1\n     -1  1]\nλ = @variable(model, lower_bound=0)\n# Symmetric matrix of nonnegative entries\nΛ = [0 λ\n     λ 0]\nusing LinearAlgebra # For `Symmetric`\n@constraint(model, Symmetric(Q - Λ) in PSDCone())\n\nwhich is equivalent to\n\n@polyvar x y\nλ = @variable(model, lower_bound=0)\n@constraint(model, x^2 - 2x*y + y^2 - 2*λ * x*y in SOSCone())\n\nwhich is the same as, using the domain keyword,\n\n@polyvar x y\n@constraint(model, x^2 - 2x*y + y^2 in SOSCone(), domain = @set x*y ≥ 0)\n\nFor consistency with its equivalent forms, the GramMatrixAttribute for this constraint is given by the gram matrix in the psd_inner cone, i.e. which should be equal to Q - Λ.\n\n[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.\n\n[MK87] K. G. Murty and S. N. Kabadi. Some NP-complete problems in quadratic and nonlinear programming. Mathematical programming, 39:117–129, 1987.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.GramMatrix",
    "page": "Constraints",
    "title": "SumOfSquares.GramMatrix",
    "category": "type",
    "text": "struct GramMatrix{T, MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} <: MP.APL{T}\n    Q::SymMatrix{T}\n    x::MVT\nend\n\nGram matrix x^top Q x where Q is a symmetric matrix indexed by the vector of monomials x.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.GramMatrixAttribute",
    "page": "Constraints",
    "title": "SumOfSquares.GramMatrixAttribute",
    "category": "type",
    "text": "struct GramMatrixAttribute <: MOI.AbstractConstraintAttribute end\n\nA constraint attribute for the GramMatrix of a constraint, that is, the positive semidefinte matrix Q indexed by the monomials in the vector X such that X^top Q X is the sum-of-squares certificate of the constraint. The\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.gram_matrix",
    "page": "Constraints",
    "title": "SumOfSquares.gram_matrix",
    "category": "function",
    "text": "gram_matrix(cref::JuMP.ConstraintRef)\n\nReturn the GramMatrixAttribute of cref.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.MomentMatrixAttribute",
    "page": "Constraints",
    "title": "SumOfSquares.MomentMatrixAttribute",
    "category": "type",
    "text": "struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute end\n\nA constraint attribute fot the MomentMatrix of a constraint.\n\n\n\n\n\n"
},

{
    "location": "constraints/#MultivariateMoments.moment_matrix",
    "page": "Constraints",
    "title": "MultivariateMoments.moment_matrix",
    "category": "function",
    "text": "moment_matrix(cref::JuMP.ConstraintRef)\n\nReturn the MomentMatrixAttribute of cref.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.CertificateMonomials",
    "page": "Constraints",
    "title": "SumOfSquares.CertificateMonomials",
    "category": "type",
    "text": "struct CertificateMonomials <: MOI.AbstractConstraintAttribute end\n\nA constraint attribute for the monomials indexing the GramMatrixAttribute and MomentMatrixAttribute certificates.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.certificate_monomials",
    "page": "Constraints",
    "title": "SumOfSquares.certificate_monomials",
    "category": "function",
    "text": "certificate_monomials(cref::JuMP.ConstraintRef)\n\nReturn the CertificateMonomials of cref.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.LagrangianMultipliers",
    "page": "Constraints",
    "title": "SumOfSquares.LagrangianMultipliers",
    "category": "type",
    "text": "struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute end\n\nA constraint attribute fot the LagrangianMultipliers assiciated to the inequalities of the domain of a constraint. There is one multiplier per inequality and no multiplier for equalities as the equalities are handled by reducing the polynomials over the ideal they generate instead of explicitely creating multipliers.\n\n\n\n\n\n"
},

{
    "location": "constraints/#SumOfSquares.lagrangian_multipliers",
    "page": "Constraints",
    "title": "SumOfSquares.lagrangian_multipliers",
    "category": "function",
    "text": "lagrangian_multipliers(cref::JuMP.ConstraintRef)\n\nReturn the LagrangianMultipliers of cref.\n\n\n\n\n\n"
},

{
    "location": "constraints/#Reference-1",
    "page": "Constraints",
    "title": "Reference",
    "category": "section",
    "text": "Inner approximations of the PSD cone that do not require semidefinite programming:SumOfSquares.DiagonallyDominantConeTriangle\nSumOfSquares.ScaledDiagonallyDominantConeTriangleApproximations of the cone of nonnegative polynomials:SumOfSquares.NonnegPolyInnerCone\nSumOfSquares.SOSCone\nSumOfSquares.SDSOSCone\nSumOfSquares.DSOSConeApproximations of the cone of positive semidefinite polynomial matrices:SumOfSquares.PSDMatrixInnerCone\nSumOfSquares.SOSMatrixConeApproximations of the cone of convex polynomials:SumOfSquares.ConvexPolyInnerCone\nSumOfSquares.SOSConvexConeApproximations of the cone of copositive matrices:SumOfSquares.CopositiveInnerAttributesGramMatrix\nSumOfSquares.GramMatrixAttribute\ngram_matrix\nSumOfSquares.MomentMatrixAttribute\nmoment_matrix\nSumOfSquares.CertificateMonomials\ncertificate_monomials\nSumOfSquares.LagrangianMultipliers\nlagrangian_multipliers"
},

]}
