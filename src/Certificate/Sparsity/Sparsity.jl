module Sparse

import MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases
using SemialgebraicSets

include("ChordalExtensionGraph.jl")
using .ChordalExtensionGraph: ChordalCompletion, ClusterCompletion
export ChordalCompletion, ClusterCompletion

import SumOfSquares

include("sparse_putinar.jl")

end
