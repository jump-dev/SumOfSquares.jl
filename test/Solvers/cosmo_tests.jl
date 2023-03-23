include("solver_preamble.jl")
using COSMO
# sosdemo9: 33000 < max_iter < 34000
factory = optimizer_with_attributes(COSMO.Optimizer, "verbose" => false, "max_iter" => 34000)
config = MOI.Test.Config(atol = 1e-3, rtol = 1e-3)
@testset "Linear" begin
    Tests.linear_test(factory, config, [
        "dsos_options_pricing",
    ])
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        "sdsos_options_pricing",
    ])
end
@testset "SDP" begin
    Tests.sd_test(
        factory,
        config,
        [
            "sos_options_pricing",
            # Iteration limit
            "chebyshev",
            # Expression: JuMP.termination_status(model) == MOI.INFEASIBLE
            # Evaluated: MathOptInterface.OPTIMAL == MathOptInterface.INFEASIBLE
            "maxcut",
        ],
    )
end
