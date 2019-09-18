module Constraint

using LinearAlgebra

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import MultivariatePolynomials
const MP = MultivariatePolynomials
import SemialgebraicSets
import MultivariateMoments
import PolyJuMP
import SumOfSquares
const SOS = SumOfSquares
const Certificate = SOS.Certificate

# Symmetric PSD matrix bridges
include("empty.jl")
include("psd2x2.jl")
include("diagonally_dominant.jl")
include("scaled_diagonally_dominant.jl")

# SOS polynomial bridges
include("utilities.jl")
include("sos_polynomial.jl")
include("sos_polynomial_in_semialgebraic_set.jl")

end
