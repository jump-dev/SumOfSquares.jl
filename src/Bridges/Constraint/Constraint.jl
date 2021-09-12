module Constraint

using LinearAlgebra

using MutableArithmetics
const MA = MutableArithmetics

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
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


# TODO bridges should redirect to `MOI.get_fallback` as well so that
# we can just use `Union{MOI.ConstraintIndex,MOI.Bridges.AbstractBridge}` in the `get_fallback` in `attributes.jl`
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.SOSDecompositionAttribute,
    bridge::Union{SOSPolynomialBridge,SOSPolynomialInSemialgebraicSetBridge},
)
    gram = MOI.get(model, SOS.GramMatrixAttribute(attr.result_index), bridge)
    return SOS.SOSDecomposition(gram, attr.ranktol, attr.dec)
end

end
