include("solver_preamble.jl")
import SDPNAL
factory = optimizer_with_attributes(
    SDPNAL.Optimizer,
    "printlevel" => 0,
    "tol" => 1e-4,
)
config = MOI.Test.Config(atol = 5e-3, rtol = 5e-3)
@testset "Linear" begin
    Tests.linear_test(factory, config, [
        # Expression: ≈(JuMP.objective_value(model), expected, atol=atol, rtol=rtol)
        # Evaluated: 142.57342982938655 ≈ 132.63 (atol=0.05, rtol=0.05)
        "dsos_options_pricing",
        # Infeasible not supported
        "dsos_horn",
    ])
end
@testset "SDP" begin
    Tests.sd_test(
        factory,
        config,
        [
            # MathOptInterface.NUMERICAL_ERROR == MathOptInterface.OPTIMAL
            "chebyshev",
            # SDPNAL supports LessThan on Nonnegatives variables but cannot
            # create free variables. See
            # https://github.com/jump-dev/MathOptInterface.jl/issues/987
            "sosdemo5_feasible",
            # Expression: ≈(JuMP.objective_value(model), expected, atol=atol, rtol=rtol)
            # Evaluated: 20.161198836088243 ≈ 17.17 (atol=0.05, rtol=0.05)
            "sos_options_pricing",
            # Infeasible not supported
            "quartic_ideal_rem",
            "quartic_ideal_2_rem",
            "sos_horn",
            "maxcut",
            "choi",
            "motzkin",
            "quadratic_infeasible_lyapunov_switched_system",
            "quadratic_infeasible_scaled_lyapunov_switched_system",
            "sosdemo5_infeasible",
            "quartic_infeasible_lyapunov_switched_system",
            "quartic_infeasible_scaled_lyapunov_switched_system",
        ],
    )
end
