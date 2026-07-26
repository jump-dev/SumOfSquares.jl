# Optimizer

Optimizers provide allows solving Polynomial Optimization programs using a hierarchy of relaxations.

```@docs
PolyJuMP.SAGE.Optimizer
SumOfSquares.Optimizer
PolyJuMP.AbstractRelaxationOptimizer
PolyJuMP.AbstractPolynomialOptimizer
```

The degree of the hierarchy is controlled by the following attribute:

```@docs
PolyJuMP.MultiplierMaxdegree
```
