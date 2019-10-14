include("solver_preamble.jl")
import SCS
factory = with_optimizer(SCS.Optimizer, verbose=0)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        # Expression: termination_status(model) == MOI.INFEASIBLE
        # Evaluated: MathOptInterface.OPTIMAL == MathOptInterface.INFEASIBLE
        "sos_horn"
    ])
end
