module Sparsity

import MultivariatePolynomials as MP
import MultivariateBases as MB
using SemialgebraicSets

include("ChordalExtensionGraph.jl")
using .ChordalExtensionGraph: ChordalCompletion, ClusterCompletion
export ChordalCompletion, ClusterCompletion

import SumOfSquares

export SignSymmetry
abstract type Pattern end
struct NoPattern <: Pattern end

include("xor_space.jl")
include("sign.jl")
include("variable.jl")
include("monomial.jl")

include("preorder.jl")
include("ideal.jl")

end
