include("solver_preamble.jl")
import SDPNAL
factory = with_optimizer(SDPNAL.Optimizer, printlevel=0, tol=1e-4)
config = MOI.Test.TestConfig(atol=5e-3, rtol=5e-3, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config, [
        # Infeasible not supported
        "dsos_horn"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        # SDPNAL supports LessThan on Nonnegatives variables but cannot
        # create free variables. See
        # https://github.com/JuliaOpt/MathOptInterface.jl/issues/987
        "sosdemo5_feasible",
        # Infeasible not supported
        "sos_horn",
        "maxcut",
        "choi_term",
        "motzkin",
        "quadratic_infeasible_lyapunov_switched_system",
        "sosdemo5_infeasible",
        "quartic_infeasible_lyapunov_switched_system"
    ])
end
