var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SumOfSquares-–-Sum-of-Squares-Programming-for-Julia-1",
    "page": "Index",
    "title": "SumOfSquares –- Sum of Squares Programming for Julia",
    "category": "section",
    "text": "SumOfSquares implements Sum of Squares reformulation for PolyJuMP, enabling the creation of sum of squares variables and constraints in JuMP.The polynomial can be represented by any library implementing the MultivariatePolynomial.jl interface. That is, you can currently choose between DynamicPolynomials and TypedPolynomials. As a rule of thumb, if you know at compile time (or at the time you are writing your code), the number of variable and that this number is small, use TypedPolynomials, otherwise, use DynamicPolynomials.Some presentations on, or using, SumOfSquares:Benoît Legat at the JuMP Meetup 2017 [Slides] [Video]\nJoey Huchette at SIAM Opt 2017The following example shows how to find lower bounds for the Goldstein-Price function using this package.using MultivariatePolynomials\nusing JuMP\nusing PolyJuMP\nusing SumOfSquares\nusing DynamicPolynomials\nusing Mosek\n\n# Create symbolic variables (not JuMP decision variables)\n@polyvar x1 x2\n\n# Create a Sum of Squares JuMP model with the Mosek solver\nm = SOSModel(solver = MosekSolver())\n\n# Create a JuMP decision variable for the lower bound\n@variable m γ\n\n# f(x) is the Goldstein-Price function\nf1 = x1+x2+1\nf2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2\nf3 = 2*x1-3*x2\nf4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2\n\nf = (1+f1^2*f2)*(30+f3^2*f4)\n\n# Constraints f(x) - γ to be sum of squares\n@constraint m f >= γ\n\n@objective m Max γ\n\nstatus = solve(m)\n\n# The lower bound found is 3\nprintln(getobjectivevalue(m))"
},

{
    "location": "index.html#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"sumofsquares.md\", \"variables.md\", \"constraints.md\"]\nDepth = 2"
},

{
    "location": "sumofsquares.html#",
    "page": "Sum-of-Squares Programming",
    "title": "Sum-of-Squares Programming",
    "category": "page",
    "text": ""
},

{
    "location": "sumofsquares.html#Sum-of-Squares-Programming-1",
    "page": "Sum-of-Squares Programming",
    "title": "Sum-of-Squares Programming",
    "category": "section",
    "text": "This section contains a brief introduction to Sum-of-Squares Programming. For more details, see [BPT12, Las09, Lau09]."
},

{
    "location": "sumofsquares.html#Quadratic-forms-and-Semidefinite-programming-1",
    "page": "Sum-of-Squares Programming",
    "title": "Quadratic forms and Semidefinite programming",
    "category": "section",
    "text": "The positive semidefiniteness of a n times n real symmetric matrix Q is equivalent to the nonnegativity of the quadratic form p(x) = x^top Q x for all vector x in mathbbR^n. For instance, the polynomialx_1^2 + 2x_1x_2 + 5x_2^2 + 4x_2x_3 + x_3^2 = x^top beginpmatrix1  1  01  5  2 0  2  1endpmatrix xis nonnegative since the matrix of the right-hand side is positive semidefinite. Moreover, a certificate of nonnegativity can be extracted from the Cholesky decomposition of the matrix:(x_1 + x_2)^2 + (2x_2 + x_3)^2 = x^top beginpmatrix1  1  00  2  1endpmatrix^top beginpmatrix1  1  00  2  1endpmatrix x"
},

{
    "location": "sumofsquares.html#Polynomial-nonnegativity-and-Semidefinite-programming-1",
    "page": "Sum-of-Squares Programming",
    "title": "Polynomial nonnegativity and Semidefinite programming",
    "category": "section",
    "text": "This can be generalized to a polynomial of arbitrary degree. A polynomial p(x) is nonnegative is it can be rewritten as p(x) = X^top Q X where Q is a real symmetric positive semidefinite matrix and X is a vector of monomials.For instancex_1^2 + 2x_1^2x_2 + 5x_1^2x_2^2 + 4x_1x_2^2 + x_2^2 = X^top beginpmatrix1  1  01  5  2 0  2  1endpmatrix Xwhere X = (x_1 x_1x_2 x_2) Similarly to the previous section, the Cholesky factorization of the matrix can be used to extract a sum of squares certificate of nonnegativity for the polynomial:(x_1 + x_1x_2)^2 + (2x_1x_2 + x_2)^2 = X^top beginpmatrix1  1  00  2  1endpmatrix^top beginpmatrix1  1  00  2  1endpmatrix X"
},

{
    "location": "sumofsquares.html#When-is-nonnegativity-equivalent-to-sum-of-squares-?-1",
    "page": "Sum-of-Squares Programming",
    "title": "When is nonnegativity equivalent to sum of squares ?",
    "category": "section",
    "text": "Determining whether a polynomial is nonnegative is NP-hard. The condition of the previous section was only sufficient, that is, there exists nonnegative polynomials that are not sums of squares. Hilbert showed in 1888 that there are exactly 3 cases for which there is equivalence between the nonnegativity of the polynomials of n variables and degree 2d and the existence of a sum of squares decomposition.n = 1 : Univariate polynomials\n2d = 2 : Quadratic polynomials\nn = 2, 2d = 4 : Bivariate quarticsThe first explicit example of polynomial that was not a sum of squares was given by Motzkin in 1967:x_1^4x_2^2 + x_1^2x_2^4 + 1 - 3x_1^2x_2^2 geq 0 quad forall xWhile it is not a sum of squares, it can still be certified to be nonnegative using sum-of-squares programming by identifying it with a rational sum-of-squares decomposition. These facts can be verified numerically using this package as detailed in the motzkin notebook."
},

{
    "location": "sumofsquares.html#References-1",
    "page": "Sum-of-Squares Programming",
    "title": "References",
    "category": "section",
    "text": "[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.[Las09] Lasserre, J. B. Moments, positive polynomials and their applications World Scientific, 2009[Lau09] Laurent, M. Sums of squares, moment matrices and optimization over polynomials Emerging applications of algebraic geometry, Springer, 2009, 157-270"
},

{
    "location": "variables.html#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": ""
},

{
    "location": "variables.html#Polynomial-and-Sum-of-Squares-variables-1",
    "page": "Variables",
    "title": "Polynomial and Sum-of-Squares variables",
    "category": "section",
    "text": ""
},

{
    "location": "variables.html#Polynomial-variables-1",
    "page": "Variables",
    "title": "Polynomial variables",
    "category": "section",
    "text": "While JuMP allows to create decision variables representing a number whose value needs to be optimized upon by the optimizer, PolyJuMP allows to create polynomial decision variables. In order to do that, we first need to create polynomial variables with the @polyvar macro:using DynamicPolynomials # or TypedPolynomials, pick your favorite\n@polyvar x yNote that these should not be confused with JuMP\'s decision variables which are created using the @variable macro. The polynomial decision variable that will be created need to be parametrized by a polynomial basis of finite size. For instance, if we want to find a quadratic polynomial, we can parametrize it with all monomials of degree between 0 and 2. Generating a vector with such monomials can be achieved through the monomials function:using MultivariatePolynomials\nX = monomials([x, y], 0:2)We can now create our polynomial variable p as follows:using PolyJuMP\n@variable(model, p, Poly(X))This creates a vector of decision variables a and sets p as the scalar product between a and X.Just like with classical JuMP\'s decision variables, containers of polynomial variables can be created as follows:using PolyJuMP\n@variable(model, p[1:3, 1:4], Poly(X))       # Creates a Matrix\n@variable(model, p[[:a, :b], -2:2], Poly(X)) # Creates a JuMPArray\n@variable(model, p[i=1:3, j=i:3], Poly(X))   # Creates a DictFor more flexibility, polynomials parametrized by decision variables can also be created \"by hand\" for instance as follows:@variable(model, α)\n@variable(model, β)\np = α*x^2 + (α+β)*y^2*x + β*y^3The coefficients of the polynomial can even be quadratic, e.g.@variable(model, α)\n@variable(model, β)\np = (3α^2+β)*x^2 + (α*β+2β)*y^2*x + β*y^3Let me stress again the distinction between α and β which are decision variables and x and y which are polynomial variables."
},

{
    "location": "variables.html#Nonnegative-polynomial-variables-1",
    "page": "Variables",
    "title": "Nonnegative polynomial variables",
    "category": "section",
    "text": "In order to create a sum-of-squares polynomial variable, the syntax is exactly the same except SOSPoly should be used instead of Poly. For instance, the following code creates a 3 times 4 matrix of sum-of-squares polynomial variables:using SumOfSquares\n@variable(model, p[1:3, 1:4], SOSPoly(X))There is however an important difference between the signification of the vector of monomials X between Poly and SOSPoly. For SOSPoly, it creates a positive semidefinite matrix of variables Q and sets p as the value of X\' * Q * X. That is, for instance, if X contains all the monomials of degree 2, then all the monomials of p will have degree 4 (i.e. p will be a quartic form).Similarly, to create diagonally-dominant-sum-of-squares polynomial variables (see [Definition 2, AM17]), use DSOSPoly(X). This creates a diagonally dominant matrix of variables Q and sets the polynomial variables as the value of X\' * Q * X.Finally, to create scaled-diagonally-dominant-sum-of-squares polynomial variables (see [Definition 2, AM17]), use DSOSPoly(X). This creates a scaled diagonally dominant matrix of variables Q and sets the polynomial variables as the value of X\' * Q * X."
},

{
    "location": "variables.html#Choosing-a-polynomial-basis-1",
    "page": "Variables",
    "title": "Choosing a polynomial basis",
    "category": "section",
    "text": "In the previous section, we show how to create polynomial variables from a monomial basis. However, the monomial basis is only a particular case of polynomial basis and while using a different basis of the same space of polynomial is would give an equivalent program, it might be more stable numerically (see [Section 3.1.5, BPT12]).For instance, creating an univariate cubic polynomial variable p using the Chebyshev basis can be done as follows:using PolyJuMP\ncheby_basis = FixedPolynomialBasis([1, x, 2x^2-1, 4x^3-3x])\n@variable(model, p, Poly(cheby_basis))and to create a quadratic form variable q using the scaled monomial basis (see [Section 3.1.5, BPT12]), use the following:using MultivariatePolynomials\nX = monomials([x], 2)\nusing PolyJuMP\nscaled_basis = ScaledMonomialBasis(X)\n@variable(model, q, Poly(scaled_basis))which is equivalent tousing PolyJuMP\nscaled_basis = FixedPolynomialBasis([x^2, √2*x*y, y^2])\n@variable(model, q, Poly(scaled_basis))"
},

{
    "location": "variables.html#References-1",
    "page": "Variables",
    "title": "References",
    "category": "section",
    "text": "[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017"
},

{
    "location": "constraints.html#",
    "page": "Constraints",
    "title": "Constraints",
    "category": "page",
    "text": ""
},

{
    "location": "constraints.html#Constraints-1",
    "page": "Constraints",
    "title": "Constraints",
    "category": "section",
    "text": ""
},

{
    "location": "constraints.html#Equality-constraints-between-polynomials-1",
    "page": "Constraints",
    "title": "Equality constraints between polynomials",
    "category": "section",
    "text": "Equality between polynomials in PolyJuMP uses the same syntax as equality between affine or quadratic expression in JuMP. For instance, creating two quadratic n-variate polynomials p and q that must sum up to one can be done as follows:using DynamicPolynomials\n@polyvar x[1:n]\nusing MultivariatePolynomials\nX = monomials(x, 0:2)\nusing PolyJuMP\n@variable(model, p, Poly(X))\n@variable(model, q, Poly(X))\n@constraint(model, p + q == 1)Vectorized constraints can also be used as well as vector of constraints, named constraints, ... For instance instance, if P and Q are two n times n matrices of polynomials, the following constraints the sum of rows and columns to match@constraint(model, con[i=1:n], sum(P[i, :]) == sum(Q[:, i]))and con[i] contains the reference to the constraint between the ith row of P and the ith column of Q."
},

{
    "location": "constraints.html#Inequality-constraints-between-polynomials-1",
    "page": "Constraints",
    "title": "Inequality constraints between polynomials",
    "category": "section",
    "text": "Polynomials can be constrained to be sum-of-squares with the in syntax. For instance, to constraints a polynomial p to be sum-of-squares, do@constraint(model, p in SOSCone())"
},

{
    "location": "constraints.html#Automatically-interpreting-polynomial-nonnegativity-as-a-sum-of-squares-constraint-1",
    "page": "Constraints",
    "title": "Automatically interpreting polynomial nonnegativity as a sum-of-squares constraint",
    "category": "section",
    "text": "As detailed in When is nonnegativity equivalent to sum of squares ?, the nonnegativity of a polynomial is not equivalent to the existence of a sum-of-squares decomposition. However, if explicitely specified, nonnegativity constraints can be automatically interpreted as sum-of-squares constraints. The simplest way to do that is to create the model withmodel = SOSModel(...)instead ofmodel = Model(...)An alternative equivalent way is to call setpolymodule! after creating the model:setpolymodule!(model, SumOfSquares)This second approach may be useful if the SumOfSquares JuMP extension need to be used with another JuMP extension that also has a special model constructor. A third alternative is the following:PolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, SOSCone)\nPolyJuMP.setdefault!(model, PolyJuMP.NonNegPolyMatrix, SOSMatrixCone)This approach adds the flexibility to choose the default cone forconstraints of the form @constraint(mode, ..., some_polynomial ≥ other_polynomial, ...) which is the cone given as default to PolyJuMP.NonNegPoly; and\nconstraints of the form @constraint(mode, ..., some_matrix_of_polynomial in PSDCone(), ...) or @SDconstraint(mode, ..., some_matrix_of_polynomial ⪰ other_matrix_of_polynomial, ...) which is the cone given as default to PolyJuMP.NonNegPolyMatrix.For instance, to use the diagonally-dominant-sum-of-squares cone (see [Definition 2, AM17]) for the first type of contraints, doPolyJuMP.setdefault!(model, PolyJuMP.NonNegPoly, DSOSCone)"
},

{
    "location": "constraints.html#Changing-the-polynomial-basis-1",
    "page": "Constraints",
    "title": "Changing the polynomial basis",
    "category": "section",
    "text": "As introduced in Choosing a polynomial basis, there may be numerical advantages to use another basis than the standard monomial basis when creating polynomial variables. Similarly, other polynomial basis can be used for polynomial constraints. However, for constraints, the polynomial space is determined by the polynomial constrained to be nonnegative. For instance, consider the constraint:@constraint(model, α * x^2 + β * y^2 ≥ (α - β) * x * y)where α and β are JuMP decision variables and x and y are polynomial variables. Since the polynomial is a quadratic form, the sum-of-squares certificate is also a quadratic form (see [Section~3.3.4, BPT12]). Hence the default polynomial basis used for the [Nonnegative polynomial variables] certificate is MonomialBasis([x, y]), that is, we search for a positive semidefinite matrix Q such that * x^2 +  * y^2 - ( - ) * x * y = X^top Q Xwhere X = (x y).As the polynomial space is determined by the polynomial being constrained, only the basis type need to be given.For instance, to use the scaled monomial basis, use@constraint(model, p ≥ q, basis = ScaledMonomialBasis)"
},

{
    "location": "constraints.html#Polynomial-nonnegativity-on-a-subset-of-the-space-1",
    "page": "Constraints",
    "title": "Polynomial nonnegativity on a subset of the space",
    "category": "section",
    "text": "By default, the constraint@constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α)constrains the polynomial to be nonnegative for every real numbers x and y. However, the set of points (x, y) for which the polynomial is constrained to be nonnegative can be specified by the domain keyword:using SemialgebraicSets\nS = @set x >= 0 && y >= 0 && x + y >= 1\n@constraint(model, x^3 - x^2 + 2x*y -y^2 + y^3 >= α, domain = S)See this notebook for a detailed example."
},

{
    "location": "constraints.html#Dual-of-polynomial-constraints-1",
    "page": "Constraints",
    "title": "Dual of polynomial constraints",
    "category": "section",
    "text": "The dual of a polynomial constraint cref is a moment serie μ as defined in MultivariateMoments. The dual can be obtained with the JuMP.resultdual function as with classical dual values in JuMP. The matrix of moment can be obtaine as follows:μ = JuMP.resultdual(cref)\nν = matmeasure(μ, certificate_monomials(cref))The extractatoms function of MultivariateMoments can be used to check if there exists an atomic measure (i.e. a measure that is a sum of Dirac measures) for that has the moments given in ν. This can be used for instance in polynomial optimization (see this notebook) or stability analysis (see this notebook)."
},

{
    "location": "constraints.html#References-1",
    "page": "Constraints",
    "title": "References",
    "category": "section",
    "text": "[BPT12] Blekherman, G.; Parrilo, P. A. & Thomas, R. R. Semidefinite Optimization and Convex Algebraic Geometry. Society for Industrial and Applied Mathematics, 2012.[AM17] Ahmadi, A. A. & Majumdar, A. DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization ArXiv e-prints, 2017"
},

]}
