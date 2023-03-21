include("solver_preamble.jl")
import GLPK
factory = GLPK.Optimizer
config = MOI.Test.Config(atol = 1e-5, rtol = 1e-5)
@testset "Linear" begin
    # With `dsos_horn_test`, the termination status is
    # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
    Tests.linear_test(factory, config, ["dsos_horn"])
end
