include("solver_preamble.jl")
import SeDuMi
factory = with_optimizer(SeDuMi.Optimizer, fid=0)
config = MOI.Test.Config(atol=1e-4, rtol=1e-4)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        # Expression: termination_status(model) == MOI.INFEASIBLE
        # Evaluated: MathOptInterface.NUMERICAL_ERROR == MathOptInterface.INFEASIBLE
        "sdsos_horn"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        "quadratic_infeasible_lyapunov_switched_system",
        "quadratic_feasible_lyapunov_switched_system",
        # ALMOST_OPTIMAL
        "maxcut", "chebyshev", "sos_quartic_comparison",
        #   Expression: JuMP.primal_status(model) == MOI.FEASIBLE_POINT
        #    Evaluated: MathOptInterface.NEARLY_FEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
        "sosdemo5_feasible"
    ])
end
