using JuMP
using SumOfSquares
using FactCheck

#x = Var("x")
#y = Var("y")
#@polyvar x y z

#println([x^2, x*y, y^2])
#println([x, y])

#println(MonomialVector((x*y*z).vars, 3))

include("solvers.jl")

include("sosdemo1.jl")
include("motzkin.jl")
include("sosdemo2.jl")
#include("sosdemo3.jl")
include("sosdemo4.jl")

FactCheck.exitstatus()
