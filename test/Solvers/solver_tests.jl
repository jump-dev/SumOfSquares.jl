import Pkg
using Test

function solver_test(name::Symbol)
    if any(dep -> dep.name == string(name), values(Pkg.dependencies()))
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
        @warn(
            "The solver $name is not part of the dependencies, run `] add $name`."
        )
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
solver_test(:SDPAFamily)
solver_test(:SDPNAL)

# SDP+SOC solvers
solver_test(:CDCS)
solver_test(:COSMO)
solver_test(:Mosek)
solver_test(:ProxSDP)
solver_test(:SeDuMi)
solver_test(:SCS)

# If we re-run `solver_tests.jl`, it may be because we changed `Tests/Tests.jl`.
shared_preamble = false
nothing # Show nothing when `include` is called from REPL
