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

]}
