using JuMP
using SumOfSquares
using FactCheck

include("solvers.jl")

include("sosdemo1.jl")
include("motzkin.jl")
include("sosdemo2.jl")
include("sosdemo3.jl")
include("sosdemo4.jl")

FactCheck.exitstatus()
