using JuMP
using SumOfSquares
using FactCheck

include("solvers.jl")

include("motzkin.jl")
include("choi.jl")

# SOSTools demos
include("sosdemo1.jl")
include("sosdemo2.jl")
include("sosdemo3.jl")
include("sosdemo4.jl")
#include("sosdemo5.jl")
include("sosdemo6.jl")
include("sosdemo9.jl")
include("sosdemo10.jl")

FactCheck.exitstatus()
