using JuMP
using SumOfSquares
using FactCheck

include("solvers.jl")

#include("motzkin.jl")

# SOSTools demos
#include("sosdemo1.jl")
include("sosdemo2.jl")
#include("sosdemo3.jl")
#include("sosdemo4.jl")
#include("sosdemo5.jl")

FactCheck.exitstatus()
