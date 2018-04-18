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
    "text": "SumOfSquares implements Sum of Squares reformulation for PolyJuMP, enabling the creation of sum of squares variables and constraints in JuMP.The polynomial can be represented by any library implementing the MultivariatePolynomial.jl interface. That is, you can currently choose between DynamicPolynomials and TypedPolynomials. As a rule of thumb, if you know at compile time (or at the time you are writing your code), the number of variable and that this number is small, use TypedPolynomials, otherwise, use DynamicPolynomials.Some presentations on, or using, SumOfSquares:Benoit Legat at the JuMP Meetup 2017\nJoey Huchette at SIAM Opt 2017The following example shows how to find lower bounds for the Goldstein-Price function using this package.using MultivariatePolynomials\nusing JuMP\nusing PolyJuMP\nusing SumOfSquares\nusing DynamicPolynomials\nusing Mosek\n\n# Create symbolic variables (not JuMP decision variables)\n@polyvar x1 x2\n\n# Create a Sum of Squares JuMP model with the Mosek solver\nm = SOSModel(solver = MosekSolver())\n\n# Create a JuMP decision variable for the lower bound\n@variable m γ\n\n# f(x) is the Goldstein-Price function\nf1 = x1+x2+1\nf2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2\nf3 = 2*x1-3*x2\nf4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2\n\nf = (1+f1^2*f2)*(30+f3^2*f4)\n\n# Constraints f(x) - γ to be sum of squares\n@constraint m f >= γ\n\n@objective m Max γ\n\nstatus = solve(m)\n\n# The lower bound found is 3\nprintln(getobjectivevalue(m))"
},

{
    "location": "index.html#Contents-1",
    "page": "Index",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"installation.md\", \"representation.md\", \"polyhedron.md\", \"redundancy.md\", \"projection.md\", \"optimization.md\", \"utilities.md\"]\nDepth = 2"
},

]}
