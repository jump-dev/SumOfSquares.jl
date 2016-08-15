using SCS
using JuMP
using SumOfSquares
using Base.Test

#x = Var("x")
#y = Var("y")
#@polyvar x y z

#println([x^2, x*y, y^2])
#println([x, y])

#println(MonomialVector((x*y*z).vars, 3))

include("sosdemo1.jl")
include("motzkin.jl")
include("sosdemo2.jl")
