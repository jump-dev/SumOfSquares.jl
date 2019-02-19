include("Tests/Tests.jl")

using Test, JuMP, CSDP
factory = with_optimizer(CSDP.Optimizer, printlevel=0)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config)
end
