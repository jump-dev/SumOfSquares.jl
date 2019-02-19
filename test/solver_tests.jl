import Pkg

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
#solver_test(:SCS)
