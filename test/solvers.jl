# Similar to JuMP/test/solvers.jl

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end

mos = try_import(:Mosek)
scs = try_import(:SCS)

# Semidefinite solvers
sdp_solvers = Any[]
# Need Mosek 8 for sosdemo3 to pass
mos && push!(sdp_solvers, Mosek.MosekSolver(LOG=0))
# Need 54000 iterations for sosdemo3 to pass
scs && push!(sdp_solvers, SCS.SCSSolver(eps=1e-6, max_iters=55000, verbose=0))
