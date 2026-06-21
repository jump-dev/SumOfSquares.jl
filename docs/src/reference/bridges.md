# Bridges

[Bridges](https://jump.dev/MathOptInterface.jl/stable/submodules/Bridges/overview/)
are the mechanism by which a Sum-of-Squares constraint is rewritten until it
hits a set that the underlying solver supports natively. See the
[Bridging mechanism](@ref) section of the manual for a high-level walk-through
of the bridge chain.

## Utilities

Bridges are automatically added by JuMP using the following PolyJuMP utilities:
```@docs
SumOfSquares.PolyJuMP.bridgeable
SumOfSquares.PolyJuMP.bridges
```

## Bridges for polynomial optimization

```@docs
PolyJuMP.ScalarPolynomialFunction
PolyJuMP.Bridges.Objective.ToPolynomialBridge
PolyJuMP.Bridges.Constraint.ToPolynomialBridge
```

## Constraint bridges

The bridges that act on the constraint side of a Sum-of-Squares constraint.

```@docs
SumOfSquares.Bridges.Constraint.SOSPolynomialBridge
SumOfSquares.Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge
SumOfSquares.Bridges.Constraint.ImageBridge
SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge
SumOfSquares.Bridges.Constraint.DiagonallyDominantBridge
SumOfSquares.Bridges.Constraint.EmptyBridge
```

## Variable bridges

The bridges that act on the variable side — they back a `WeightedSOSCone` or
one of its PSD-approximation cones with constrained variables.

```@docs
SumOfSquares.Bridges.Variable.KernelBridge
SumOfSquares.Bridges.Variable.LowRankBridge
SumOfSquares.Bridges.Variable.PositiveSemidefinite2x2Bridge
SumOfSquares.Bridges.Variable.ScaledDiagonallyDominantBridge
SumOfSquares.Bridges.Variable.CopositiveInnerBridge
```

## SAGE bridges

```@docs
SumOfSquares.PolyJuMP.SAGE.SignomialsBridge
SumOfSquares.PolyJuMP.SAGE.AGEBridge
```
