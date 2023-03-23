include("solver_preamble.jl")
import CSDP
factory = optimizer_with_attributes(CSDP.Optimizer, "printlevel" => 0)
config = MOI.Test.Config(atol = 1e-4, rtol = 1e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config)
end
