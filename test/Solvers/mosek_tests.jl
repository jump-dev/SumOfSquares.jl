include("solver_preamble.jl")
using MosekTools
# It cannot be used in Direct mode yet as it gets:
# MathOptInterface.DeleteNotAllowed{MathOptInterface.ConstraintIndex{MathOptInterface.VectorOfVariables,MathOptInterface.RotatedSecondOrderCone}}
# for `sdsos_bivariate_quadratic`
factory = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
config = MOI.Test.Config(atol = 1e-5, rtol = 1e-5)
@testset "Linear" begin
    Tests.linear_test(factory, config)
end
@testset "SOC" begin
    Tests.soc_test(factory, config)
end
@testset "SDP" begin
    Tests.sd_test(
        factory,
        config,
        [
            # Expression: termination_status(model) == MOI.OPTIMAL
            #  Evaluated: MathOptInterface.SLOW_PROGRESS == MathOptInterface.OPTIMAL
            # Expression: ≈(objective_value(model), α_value, atol=atol, rtol=rtol)
            #  Evaluated: 0.6847213936286884 ≈ 10.0 (atol=1.0e-5, rtol=1.0e-5)
            "BPT12e399_maxdegree",
            # FIXME AssertionError: m.x_type[ref2id(vi)] == Deleted
            "maxcut",
            "motzkin",
            # Expression: JuMP.termination_status(model) == MOI.INFEASIBLE
            # Evaluated: MathOptInterface.SLOW_PROGRESS == MathOptInterface.INFEASIBLE
            # Contacted Mosek and they replied that there is nothing wrong as the
            # PrimalStatus and DualStatus are correct
            "quartic_infeasible_lyapunov_switched_system",
            "quartic_infeasible_scaled_lyapunov_switched_system",
        ],
    )
end
