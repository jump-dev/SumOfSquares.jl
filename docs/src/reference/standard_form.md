# Standard form

## JuMP variables

```@docs
SumOfSquares.PolyJuMP.Poly
```

## JuMP Sets

Approximations of the cone of nonnegative polynomials:
```@docs
SumOfSquares.NonnegPolyInnerCone
SumOfSquares.SOSCone
SumOfSquares.SDSOSCone
SumOfSquares.DSOSCone
```

Approximations of the cone of positive semidefinite polynomial matrices:
```@docs
SumOfSquares.PSDMatrixInnerCone
SumOfSquares.SOSMatrixCone
```

Approximations of the cone of convex polynomials:
```@docs
SumOfSquares.ConvexPolyInnerCone
SumOfSquares.SOSConvexCone
```

Approximations of the cone of copositive matrices:
```@docs
SumOfSquares.CopositiveInner
```

SAGE cones:
```@docs
SumOfSquares.PolyJuMP.SAGE.Polynomials
SumOfSquares.PolyJuMP.SAGE.Signomials
SumOfSquares.PolyJuMP.SAGE.SignomialsBridge
SumOfSquares.PolyJuMP.SAGE.AGEBridge
```

## MOI Sets

Special cases of positive semidefinite cones:
```@docs
SumOfSquares.EmptyCone
SumOfSquares.PositiveSemidefinite2x2ConeTriangle
```

Inner approximations of the PSD cone that do not require semidefinite
programming:
```@docs
SumOfSquares.DiagonallyDominantConeTriangle
SumOfSquares.ScaledDiagonallyDominantConeTriangle
SumOfSquares.Bridges.Variable.ScaledDiagonallyDominantBridge
```

## Bridges

Bridges are automatically added using the following utilities:
```@docs
SumOfSquares.PolyJuMP.bridgeable
SumOfSquares.PolyJuMP.bridges
```

Bridges for polynomial optimization
```@docs
PolyJuMP.ScalarPolynomialFunction
PolyJuMP.Bridges.Objective.ToPolynomialBridge
PolyJuMP.Bridges.Constraint.ToPolynomialBridge
```
