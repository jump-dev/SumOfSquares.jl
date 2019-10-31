include("solver_preamble.jl")
import SDPAFamily
factory = with_optimizer(SDPAFamily.Optimizer{Float64}, params=(gammaStar=0.8, maxIteration=200),presolve=true)
config = MOI.Test.TestConfig(atol=1e-5, rtol=1e-5, query=false)
@testset "Linear" begin
    # With `dsos_horn`, the termination status is
    # `INFEASIBLE_OR_UNBOUNDED` instead of `INFEASIBLE`.
    #
    # With `dsos_term`, we get the error
    # ArgumentError: Bridge of type `MathOptInterface.Bridges.Constraint.VectorFunctionizeBridge{Float64,SumOfSquares.SOSPolynomialSet{SemialgebraicSets.FullSpace,DynamicPolynomials.Monomial{true},DynamicPolynomials.MonomialVector{true},SumOfSquares.Certificate.Remainder{SumOfSquares.NonnegPolyInnerCone{SumOfSquares.DiagonallyDominantConeTriangle},PolyJuMP.MonomialBasis,Tuple{}}}}` does not support accessing the attribute `SumOfSquares.GramMatrixAttribute(1)`.
    Tests.linear_test(factory, config, ["dsos_term", "dsos_horn"])
end
@testset "SDP" begin
    # With `maxcut` and `quartic_infeasible_lyapunov_switched_system`, the dual status is `UNKNOWN_RESULT_STATUS` instead of `MOI.INFEASIBILITY_CERTIFICATE`
    #
    # With `sos_term`, we get
    # ArgumentError: Bridge of type `MathOptInterface.Bridges.Constraint.VectorFunctionizeBridge{Float64,SumOfSquares.SOSPolynomialSet{SemialgebraicSets.FullSpace,DynamicPolynomials.Monomial{true},DynamicPolynomials.MonomialVector{true},SumOfSquares.Certificate.Remainder{SumOfSquares.NonnegPolyInnerCone{MathOptInterface.PositiveSemidefiniteConeTriangle},PolyJuMP.MonomialBasis,Tuple{}}}}` does not support accessing the attribute `SumOfSquares.GramMatrixAttribute(1)`.
    Tests.sd_test(factory, config, ["maxcut", "quartic_infeasible_lyapunov_switched_system", "sos_term"])

end
