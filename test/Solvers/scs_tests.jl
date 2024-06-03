include("solver_preamble.jl")
import SCS
factory = optimizer_with_attributes(SCS.Optimizer, MOI.Silent() => true)
config = MOI.Test.Config(atol = 1e-3, rtol = 1e-3)
@testset "Linear" begin
    Tests.linear_test(factory, config, include=["dsos_horn",#"dsos_concave_then_convex_cubic"
    ])
end
#@testset "SOC" begin
#    Tests.soc_test(factory, config)
#end
#@testset "SDP" begin
#    Tests.sd_test(factory, config)
#end
