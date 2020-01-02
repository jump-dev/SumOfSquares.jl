include("solver_preamble.jl")
import SCS
factory = with_optimizer(SCS.Optimizer, verbose=0)
config = MOI.Test.TestConfig(atol=5e-4, rtol=5e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        "sdsos_options_pricing"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        "sos_options_pricing",
        # Expression: termination_status(model) == MOI.INFEASIBLE
        # Evaluated: MathOptInterface.OPTIMAL == MathOptInterface.INFEASIBLE
        "sos_horn"
    ])
end
