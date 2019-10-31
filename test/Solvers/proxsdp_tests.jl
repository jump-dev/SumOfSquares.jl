include("solver_preamble.jl")
import ProxSDP
factory = with_optimizer(ProxSDP.Optimizer, log_verbose=false)
config = MOI.Test.TestConfig(atol=2e-2, rtol=2e-2)
@testset "Linear" begin
    Tests.linear_test(factory, config, [
        # Expression: termination_status(model) == MOI.OPTIMAL
        #  Evaluated: MathOptInterface.INFEASIBLE_OR_UNBOUNDED == MathOptInterface.OPTIMAL
        "dsos_concave_then_convex_cubic",
        # Expression: JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        #  Evaluated: MathOptInterface.INFEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        "dsos_quartic_comparison",
        "dsos_horn"
    ])
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        # Expression: termination_status(model) == MOI.OPTIMAL
        #  Evaluated: MathOptInterface.INFEASIBLE_OR_UNBOUNDED == MathOptInterface.OPTIMAL
        "sdsos_concave_then_convex_cubic",
        # ConstraintDual not implemented for VoV-in-SOC
        "sdsos_univariate_quadratic",
        "sdsos_horn", "sdsos_bivariate_quadratic"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        # with γ=3.9 it should be infeasible: Test Failed at /home/blegat/.julia/dev/SumOfSquares/test/Tests/maxcut.jl:37
        # Expression: JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        #  Evaluated: MathOptInterface.INFEASIBLE_POINT == MathOptInterface.INFEASIBILITY_CERTIFICATE
        "maxcut",
        # ITERATION_LIMIT
        "choi_term", "motzkin",
        # ArgumentError: ModelLike of type ProxSDP.Optimizer does not support accessing the attribute MathOptInterface.ConstraintDual(1)
        "sos_univariate_quadratic",
        "sos_horn",
        "sos_bivariate_quadratic",
        "quadratic_infeasible_lyapunov_switched_system",
        "quartic_infeasible_lyapunov_switched_system"
    ])
end
