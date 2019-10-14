include("solver_preamble.jl")
using COSMO
# sosdemo9: 33000 < max_iter < 34000
factory = with_optimizer(COSMO.Optimizer, verbose=false, max_iter = 34000)
config = MOI.Test.TestConfig(atol=1e-3, rtol=1e-3)
@testset "Linear" begin
    Tests.linear_test(factory, config, [
        # See https://github.com/oxfordcontrol/COSMO.jl/issues/96
        "dsos_horn"
    ])
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        # See https://github.com/oxfordcontrol/COSMO.jl/issues/96
        "sdsos_horn"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        # Iteration limit
        "chebyshev",
        # Expression: JuMP.termination_status(model) == MOI.INFEASIBLE
        # Evaluated: MathOptInterface.OPTIMAL == MathOptInterface.INFEASIBLE
        "quartic_infeasible_lyapunov_switched_system",
        "maxcut", "sos_horn", "motzkin"
    ])
end
