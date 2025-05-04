module Constraint

using LinearAlgebra

import MutableArithmetics as MA
import StarAlgebras as SA

import MathOptInterface as MOI

import MultivariatePolynomials as MP
import MultivariateBases as MB
import SemialgebraicSets
import MultivariateMoments
import PolyJuMP
import SumOfSquares as SOS
const Certificate = SOS.Certificate

# Symmetric PSD matrix bridges
include("empty.jl")
include("psd2x2.jl")
include("diagonally_dominant.jl")

# SOS polynomial bridges
include("utilities.jl")
include("image.jl")
include("sos_polynomial.jl")
include("sos_polynomial_in_semialgebraic_set.jl")

function add_all_bridges(model, ::Type{T}) where {T}
    MOI.Bridges.add_bridge(model, EmptyBridge{T})
    MOI.Bridges.add_bridge(model, PositiveSemidefinite2x2Bridge{T})
    MOI.Bridges.add_bridge(model, DiagonallyDominantBridge{T})
    MOI.Bridges.add_bridge(model, ImageBridge{T})
    MOI.Bridges.add_bridge(model, SOSPolynomialBridge{T})
    MOI.Bridges.add_bridge(model, SOSPolynomialInSemialgebraicSetBridge{T})
end

end
