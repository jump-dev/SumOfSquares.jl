include("solver_preamble.jl")
import ECOS
factory = with_optimizer(ECOS.Optimizer, verbose=false)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
