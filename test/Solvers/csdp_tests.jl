include("solver_preamble.jl")
import CSDP
factory = with_optimizer(CSDP.Optimizer, printlevel=0)
config = MOI.Test.TestConfig(atol=1e-4, rtol=1e-4, query=false)
Tests.BPT12e399_rem_test(factory, config)
Tests.BPT12e399_maxdegree_test(factory, config)
#Tests.quartic_ideal_test(factory, config)
#Tests.quartic_ideal_4_test(factory, config)
#Tests.quartic_ideal_rem_test(factory, config)
#Tests.quartic_ideal_2_rem_test(factory, config)
#Tests.quartic_ideal_4_rem_test(factory, config)
#@testset "Linear" begin
#    Tests.linear_test(factory, config)
#end
#@testset "SDP" begin
#    Tests.sd_test(factory, config)
#end
