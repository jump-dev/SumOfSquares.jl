include("solver_preamble.jl")
import SDPA
factory = with_optimizer(SDPA.Optimizer)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config, [
        # The termination status is
        # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
        "dsos_horn",
        # cholesky miss condition :: not positive definite :: info = 32 :: line 785 in sdpa_linear.cpp
        # There are some possibilities. :: line 786 in sdpa_linear.cpp
        # 1. SDPA terminates due to inaccuracy of numerical error :: line 787 in sdpa_linear.cpp
        # 2. The input problem may not have (any) interior-points :: line 788 in sdpa_linear.cpp
        # 3. Input matrices are linearly dependent :: line 789 in sdpa_linear.cpp
        "dsos_concave_then_convex_cubic",
        # The termination status is
        # `SLOW_PROGRESS` instead of `INFEASIBLE`.
        # cholesky miss condition :: not positive definite :: info = 7 :: line 785 in sdpa_linear.cpp
        # There are some possibilities. :: line 786 in sdpa_linear.cpp
        # 1. SDPA terminates due to inaccuracy of numerical error :: line 787 in sdpa_linear.cpp
        # 2. The input problem may not have (any) interior-points :: line 788 in sdpa_linear.cpp
        # 3. Input matrices are linearly dependent :: line 789 in sdpa_linear.cpp
        "dsos_bivariate_quadratic", "dsos_univariate_quadratic"])
end
@testset "SDP" begin
    # With `sos_concave_then_convex_cubic`, we have the same error than with
    # `dsos_concave_then_convex_cubic`.
    Tests.sd_test(factory, config, [
        "sos_concave_then_convex_cubic",
        "quartic_infeasible_lyapunov_switched_system",
        # Strange behavior : primal < dual :: line 144 in sdpa_solve.cpp
        # pdINF criteria :: line 1192 in sdpa_parts.cpp
        "maxcut",
        # Strange behavior : primal < dual :: line 144 in sdpa_solve.cpp
        # dUNBD criteria :: line 1220 in sdpa_parts.cpp
        # dUNBD criteria :: line 1220 in sdpa_parts.cpp
        # dUNBD criteria :: line 1220 in sdpa_parts.cpp
        # Step length is too small.  :: line 198 in sdpa_dataset.cpp
        # cannot move: step length is too short :: line 176 in sdpa_solve.cpp
        "sos_horn"
    ])
end
