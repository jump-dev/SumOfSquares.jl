using JuMP
using SumOfSquares
using FactCheck

include("solvers.jl")

#include("mono.jl")
#include("poly.jl")
#include("rational.jl")
include("promote.jl")
include("comp.jl")
include("alg.jl")
include("diff.jl")
include("subs.jl")
#include("show.jl")
include("certificate.jl")
#include("jump.jl")

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
