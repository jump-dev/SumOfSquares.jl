include("solver_preamble.jl")
using MosekTools
# It cannot be used in Direct mode yet as it gets:
# MathOptInterface.DeleteNotAllowed{MathOptInterface.ConstraintIndex{MathOptInterface.VectorOfVariables,MathOptInterface.RotatedSecondOrderCone}}
# for `sdsos_bivariate_quadratic`
factory = with_optimizer(Mosek.Optimizer, QUIET=true)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config, [
        # FIXME Mosek.MosekError(1307, "Variable '' (29) is a member of cone '' (9).")
        "sdsos_horn"
    ])
end
@testset "SDP" begin
    Tests.sd_test(factory, config, [
        # FIXME AssertionError: m.x_type[ref2id(vi)] == Deleted
        "maxcut", "sos_horn", "motzkin",
        # Expression: JuMP.termination_status(model) == MOI.INFEASIBLE
        # Evaluated: MathOptInterface.SLOW_PROGRESS == MathOptInterface.INFEASIBLE
        "quartic_infeasible_lyapunov_switched_system"
    ])
end
