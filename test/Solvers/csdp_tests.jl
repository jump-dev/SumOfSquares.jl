include("solver_preamble.jl")
import CSDP
factory = optimizer_with_attributes(CSDP.Optimizer, "printlevel" => 0)
config = MOI.Test.Config(atol = 1e-4, rtol = 1e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(
        factory,
        config,
        [
            "sosdemo5_infeasible",
            "maxcut",
            "quadratic_infeasible_lyapunov_switched_system",
            "quadratic_infeasible_scaled_lyapunov_switched_system",
            "quartic_infeasible_lyapunov_switched_system",
            "motzkin",
            "quartic_infeasible_scaled_lyapunov_switched_system",
            "choi",
        ]
    )
end
