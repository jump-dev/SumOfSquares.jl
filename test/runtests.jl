using MultivariatePolynomials
using JuMP
using SumOfSquares
using PolyJuMP
using Base.Test

#include("variable.jl")
#include("constraint.jl")

#include("solvers.jl")

include("certificate.jl")

#include("motzkin.jl")

# SOSTools demos and other examples
#include("sospoly.jl")
#include("lyapunov.jl")
#include("sosdemo3.jl")
#include("sosdemo4.jl")
#include("sosdemo5.jl")
#include("sosdemo6.jl")
#include("domain.jl")
#include("sosmatrix.jl")
#include("equalitypolyconstr.jl")
