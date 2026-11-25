# Certificate

Default choice for the `maxdegree` keyword and truncation of the equalities
of the domain whose Gröbner basis is not known:
```@docs
PolyJuMP.default_maxdegree
PolyJuMP.truncated_domain
```

Types of sparsity
```@docs
SumOfSquares.Certificate.Sparsity.Variable
SumOfSquares.Certificate.Sparsity.Monomial
SumOfSquares.Certificate.Sparsity.VariableMonomial
SumOfSquares.Certificate.Sparsity.SignSymmetry
SumOfSquares.Certificate.Sparsity.XORSpace
```

Ideal certificates:
```@docs
SumOfSquares.Certificate.MaxDegree
SumOfSquares.Certificate.FixedBasis
SumOfSquares.Certificate.Newton
SumOfSquares.Certificate.Remainder
SumOfSquares.Certificate.Sparsity.Ideal
SumOfSquares.Certificate.Symmetry.Ideal
```

Preorder certificates:
```@docs
SumOfSquares.Certificate.Putinar
SumOfSquares.Certificate.Sparsity.Preorder
```

Chordal extension:
```@docs
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.neighbors
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.is_clique
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.LabelledGraph
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.add_node!
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.add_edge!
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.add_clique!
SumOfSquares.Certificate.Sparsity.ChordalExtensionGraph.completion
```

Symmetry reduction:
```@docs
SumOfSquares.Certificate.Symmetry.InvariantBasis
SumOfSquares.Certificate.Symmetry.gram_character_degrees
SumOfSquares.Certificate.gram_weights
```
