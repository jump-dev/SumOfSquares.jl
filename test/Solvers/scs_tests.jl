include("solver_preamble.jl")
import SCS
factory = optimizer_with_attributes(SCS.Optimizer, "verbose" => 0)
config = MOI.Test.Config(atol = 5e-4, rtol = 5e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config)
end
