module Variable

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import SumOfSquares
const SOS = SumOfSquares

include("psd2x2.jl")
include("scaled_diagonally_dominant.jl")
include("copositive_inner.jl")

end
