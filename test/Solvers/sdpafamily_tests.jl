include("solver_preamble.jl")
import SDPAFamily
factory = optimizer_with_attributes(
    SDPAFamily.Optimizer{Float64},
    "params" => (gammaStar = 0.8, maxIteration = 200),
    "presolve" => true,
)
config = MOI.Test.Config(atol = 1e-5, rtol = 1e-5)
@testset "Linear" begin
    # With `dsos_horn`, the termination status is
    # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
    Tests.linear_test(factory, config, ["dsos_horn"])
end
@testset "SDP" begin
    # With `maxcut` and `quartic_infeasible_lyapunov_switched_system`, the dual status is `UNKNOWN_RESULT_STATUS` instead of `MOI.INFEASIBILITY_CERTIFICATE`
    Tests.sd_test(
        factory,
        config,
        ["maxcut", "quartic_infeasible_lyapunov_switched_system"],
    )
end
