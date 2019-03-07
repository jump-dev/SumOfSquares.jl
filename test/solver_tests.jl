import Pkg
using Test

function solver_test(name::Symbol)
    if string(name) in keys(Pkg.installed())
        ok = false
        try
            @eval import $name
            ok = true
        catch e
            @warn("The solver $name cannot be imported, run `] build $name`.")
        end
        if ok
            @testset "$name" begin
                include("$(lowercase(string(name)))_tests.jl")
            end
        end
    else
        @warn("The solver $name is not installed, run `] add $name`.")
    end
end

include("solver_preamble.jl")
shared_preamble = true

# LP solvers
solver_test(:GLPK)

# SOCP solvers
solver_test(:ECOS)

# SDP solvers (SOC is reformulated into SDP)
solver_test(:CSDP)
solver_test(:SDPA)

# SDP+SOC solvers
solver_test(:Mosek)
solver_test(:SeDuMi)
solver_test(:SCS)

# If we re-run `solver_tests.jl`, it may be because we changed `Tests/Tests.jl`.
shared_preamble = false
nothing # Show nothing when `include` is called from REPL
