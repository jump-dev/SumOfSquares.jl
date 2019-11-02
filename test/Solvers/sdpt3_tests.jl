include("solver_preamble.jl")
import SDPT3
factory = with_optimizer(SDPT3.Optimizer, printlevel = 0)
config = MOI.Test.TestConfig(atol=1e-4, rtol=1e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config, String[
        #  Expression: termination_status(model) == MOI.OPTIMAL
        #   Evaluated: MathOptInterface.SLOW_PROGRESS == MathOptInterface.OPTIMAL
        #  Expression: primal_status(model) == MOI.FEASIBLE_POINT
        #   Evaluated: MathOptInterface.NO_SOLUTION == MathOptInterface.FEASIBLE_POINT
        #  Expression: ≈(value(p), x ^ 3, atol=atol, rtol=rtol)
        #   Evaluated: 1.0173770753701414x³ + 1.5745867283645243e-11x² - 5.33707366495717e-12x + 7.1247443708677416e-12 ≈ x³ (atol=0.0001, rtol=0.0001)
        "sdsos_concave_then_convex_cubic"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, String[
        #  Expression: JuMP.termination_status(model) == MOI.OPTIMAL
        #   Evaluated: MathOptInterface.OTHER_ERROR == MathOptInterface.OPTIMAL
        #  Expression: JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        #   Evaluated: MathOptInterface.NO_SOLUTION == MathOptInterface.FEASIBLE_POINT
        "chebyshev"
    ])
end
