# Constraints

## Attributes

```@docs
SumOfSquares.PolyJuMP.MomentsAttribute
SumOfSquares.MultivariateMoments.moments(::SumOfSquares.JuMP.ConstraintRef)
GramMatrix
SumOfSquares.GramMatrixAttribute
gram_matrix
gram_operate
SumOfSquares.MomentMatrixAttribute
moment_matrix
SumOfSquares.CertificateBasis
certificate_basis
certificate_monomials
SumOfSquares.LagrangianMultipliers
lagrangian_multipliers
SOSDecomposition
SOSDecompositionWithDomain
SumOfSquares.SOSDecompositionAttribute
sos_decomposition
SumOfSquares.MultiplierIndexBoundsError
```

SAGE decomposition attribute:
```@docs
SumOfSquares.PolyJuMP.SAGE.Decomposition
SumOfSquares.PolyJuMP.SAGE.DecompositionAttribute
```

## Rounding

Heuristics to find a feasible solution from the moments computed by the
Sum-of-Squares program:
```@docs
round_solution
SumOfSquares.AbstractRounding
FirstMomentRounding
GaussianRounding
rounding_candidates
heuristic_projection
```
