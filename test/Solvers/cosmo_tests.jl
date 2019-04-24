include("solver_preamble.jl")
using COSMO
factory = with_optimizer(COSMO.Optimizer)
config = MOI.Test.TestConfig(atol=1e-4, rtol=1e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config)
end
