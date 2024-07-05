import SumOfSquares as SOS
import MultivariateBases as MB
using SemialgebraicSets
using DynamicPolynomials

function generate(test, solver, config)
    optimizer = MOI.instantiate(solver, with_bridge_type = Float64)
    # We don't use `UniversalFallback` so that the SOS bridges
    # are applied
    cached = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.Model{Float64}(),
        optimizer,
    )
    bridged = MOI.Bridges.full_bridge_optimizer(cached, Float64)
    # Second cache needed for the fallback for `ConstraintPrimal`
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        bridged,
    )
    MOI.Bridges.add_bridge(
        bridged,
        SOS.Bridges.Constraint.SOSPolynomialBridge{Float64},
    )
    MOI.Bridges.add_bridge(bridged, SOS.Bridges.Variable.KernelBridge{Float64})
    MOI.Bridges.add_bridge(bridged, SOS.Bridges.Constraint.EmptyBridge{Float64})
    MOI.Bridges.add_bridge(
        bridged,
        SOS.Bridges.Constraint.PositiveSemidefinite2x2Bridge{Float64},
    )
    MOI.Bridges.add_bridge(
        bridged,
        SOS.Bridges.Constraint.DiagonallyDominantBridge{Float64},
    )
    MOI.Bridges.add_bridge(
        bridged,
        SOS.Bridges.Constraint.ScaledDiagonallyDominantBridge{Float64},
    )
    MOI.Bridges.add_bridge(
        bridged,
        SOS.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge{Float64},
    )
    return test(model, config)
end
