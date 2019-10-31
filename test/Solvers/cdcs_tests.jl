include("solver_preamble.jl")
import CDCS
# Iterations:
# dsos_concave_then_convex_cubic : > 2000, < 3000
# chebyshev : > 12000, < 12500
factory = with_optimizer(CDCS.Optimizer, verbose=0, maxIter=12500)
# chebyshev : > 2e-3, < 3e-3
config = MOI.Test.TestConfig(atol=3e-3, rtol=3e-3, query=false)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(factory, config)
end
