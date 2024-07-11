module Variable

import SparseArrays
import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariateBases as MB
import MathOptInterface as MOI
import MultivariatePolynomials as MP
import SumOfSquares as SOS

include("psd2x2.jl")
include("scaled_diagonally_dominant.jl")
include("copositive_inner.jl")
include("kernel.jl")
include("lowrank.jl")

end
