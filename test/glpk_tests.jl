include("solver_preamble.jl")
import GLPK
factory = with_optimizer(GLPK.Optimizer)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    # With `dsos_horn_test`, the termination status is
    # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
    Tests.linear_test(factory, config, ["dsos_horn"])
end
