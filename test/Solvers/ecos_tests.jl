include("solver_preamble.jl")
import ECOS
factory = optimizer_with_attributes(ECOS.Optimizer, "verbose" => false)
config = MOI.Test.Config(atol = 1e-5, rtol = 1e-5)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        # K = 30: Test Failed at /home/blegat/.julia/dev/SumOfSquares/test/Tests/options_pricing.jl:31
        #   Expression: JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        #    Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        # K = 30: Test Failed at /home/blegat/.julia/dev/SumOfSquares/test/Tests/options_pricing.jl:32
        #   Expression: ≈(JuMP.objective_value(model), expected, atol=atol, rtol=rtol)
        #    Evaluated: 21.511440549896246 ≈ 21.51 (atol=1.0e-5, rtol=1.0e-5)
        "sdsos_options_pricing",
    ])
end
