include("solver_preamble.jl")
import SDPAFamily
factory = with_optimizer(SDPAFamily.Optimizer{Float64}, params=(gammaStar=0.8, maxIteration=200),presolve=true)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    # With `dsos_horn`, the termination status is
    # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
    Tests.linear_test(factory, config, ["dsos_horn"])
end
@testset "SDP" begin
    # With `maxcut` and `quartic_infeasible_lyapunov_switched_system`, the dual status is `UNKNOWN_RESULT_STATUS` instead of `MOI.INFEASIBILITY_CERTIFICATE`
    Tests.sd_test(factory, config, ["maxcut", "quartic_infeasible_lyapunov_switched_system"])

end
