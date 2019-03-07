include("solver_preamble.jl")
import SDPA
factory = with_optimizer(SDPA.Optimizer)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    # With `dsos_horn_test`, the termination status is
    # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
    # With `dsos_concave_then_convex_cubic`, we have
    # cholesky miss condition :: not positive definite :: info = 32 :: line 785 in sdpa_linear.cpp
    # There are some possibilities. :: line 786 in sdpa_linear.cpp
    # 1. SDPA terminates due to inaccuracy of numerical error :: line 787 in sdpa_linear.cpp
    # 2. The input problem may not have (any) interior-points :: line 788 in sdpa_linear.cpp
    # 3. Input matrices are linearly dependent :: line 789 in sdpa_linear.cpp
    Tests.linear_test(factory, config,
                      ["dsos_horn", "dsos_concave_then_convex_cubic"])
end
@testset "SDP" begin
    # With `sos_concave_then_convex_cubic`, we have the same error than with
    # `dsos_concave_then_convex_cubic`.
    Tests.sd_test(factory, config, ["sos_concave_then_convex_cubic"])
end
