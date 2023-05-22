module Constraint

using LinearAlgebra

import MutableArithmetics as MA

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
include("scaled_diagonally_dominant.jl")

# SOS polynomial bridges
include("utilities.jl")
include("image.jl")
include("sos_polynomial.jl")
include("sos_polynomial_in_semialgebraic_set.jl")

# TODO bridges should redirect to `MOI.get_fallback` as well so that
# we can just use `Union{MOI.ConstraintIndex,MOI.Bridges.AbstractBridge}` in the `get_fallback` in `attributes.jl`
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.SOSDecompositionAttribute,
    bridge::Union{GeometricBridge,SOSPolynomialBridge,SOSPolynomialInSemialgebraicSetBridge},
)
    gram = MOI.get(model, SOS.GramMatrixAttribute(attr.result_index), bridge)
    return SOS.SOSDecomposition(gram, attr.ranktol, attr.dec)
end

end
