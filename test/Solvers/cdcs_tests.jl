include("solver_preamble.jl")
import CDCS
# Iterations:
# dsos_concave_then_convex_cubic : > 2000, < 3000
factory = with_optimizer(CDCS.Optimizer, verbose=0, maxIter=3000)
config = MOI.Test.TestConfig(atol=1e-3, rtol=1e-3, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config)
end
