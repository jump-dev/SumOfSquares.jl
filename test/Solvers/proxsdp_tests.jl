include("solver_preamble.jl")
import ProxSDP
factory = optimizer_with_attributes(ProxSDP.Optimizer, "log_verbose" => false)
config = MOI.Test.Config(atol = 2e-2, rtol = 2e-2)
@testset "Linear" begin
    Tests.linear_test(
        factory,
        config,
        [
            "dsos_options_pricing",
            # Expression: JuMP.primal_status(model) == MOI.FEASIBLE_POINT
            #  Evaluated: MathOptInterface.INFEASIBLE_POINT == MathOptInterface.FEASIBLE_POINT
            "dsos_quartic_comparison",
            "dsos_cheby_bivariate_quadratic",
            "dsos_cheby_univariate_quadratic",
            "dsos_term",
            "dsos_term_fixed",
            "dsos_scaled_bivariate_quadratic",
            "dsos_scaled_univariate_quadratic",
            "dsos_bivariate_quadratic",
            "dsos_univariate_quadratic",
        ],
    )
end
@testset "SOC" begin
    Tests.soc_test(
        factory,
        config,
        [
            "sdsos_options_pricing",
            # ConstraintDual not implemented for VoV-in-SOC
            "sdsos_univariate_quadratic",
            "sdsos_bivariate_quadratic",
            "sdsos_term_fixed",
            "sdsos_term",
        ],
    )
end
@testset "SDP" begin
    Tests.sd_test(
        factory,
        config,
        [
            "sos_options_pricing",
            # ProxSDP does not return duals for nonlinear Conic constraints. Only linear constraints (equalities and inequalities) can be queried.
            "quartic_ideal",
            "quartic_ideal_4",
            "quartic_ideal_4_rem",
            # with γ=3.9 it should be infeasible: Test Failed at /home/blegat/.julia/dev/SumOfSquares/test/Tests/maxcut.jl:37
            # Expression: JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
            #  Evaluated: MathOptInterface.INFEASIBLE_POINT == MathOptInterface.INFEASIBILITY_CERTIFICATE
            "maxcut",
            # ITERATION_LIMIT
            "choi",
            "motzkin",
            # ArgumentError: ModelLike of type ProxSDP.Optimizer does not support accessing the attribute MathOptInterface.ConstraintDual(1)
            "sos_univariate_quadratic",
            "sos_bivariate_quadratic",
            "quadratic_infeasible_lyapunov_switched_system",
            "quartic_infeasible_lyapunov_switched_system",
            "sosdemo5_infeasible",
            "sos_term_fixed",
            "sos_term",
        ],
    )
end
