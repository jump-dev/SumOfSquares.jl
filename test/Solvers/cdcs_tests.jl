include("solver_preamble.jl")
import CDCS
# Iterations:
# dsos_concave_then_convex_cubic : > 2000, < 3000
# chebyshev : > 12000, < 12500
factory = optimizer_with_attributes(CDCS.Optimizer, "verbose" => 0, "maxIter" => 12500)
# chebyshev : > 2e-3, < 3e-3
config = MOI.Test.Config(atol=3e-3, rtol=3e-3)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        "sdsos_options_pricing"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        "sos_options_pricing"
    ])
end
