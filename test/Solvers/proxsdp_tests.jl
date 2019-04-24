include("solver_preamble.jl")
import ProxSDP
factory = with_optimizer(ProxSDP.Optimizer, verbose=0)
config = MOI.Test.TestConfig(atol=1e-3, rtol=1e-3)
@testset "Linear" begin
    Tests.linear_test(factory, config, ["dsos_horn"])
end
@testset "SOC" begin
    Tests.soc_test(factory, config, ["sdsos_horn", "sdsos_bivariate_quadratic"])
end
@testset "SDP" begin
    Tests.sd_test(factory, config,
                  ["sos_horn",
                   "sos_bivariate_quadratic",
                   "quadratic_infeasible_lyapunov_switched_system",
                   "quartic_infeasible_lyapunov_switched_system"])
end
