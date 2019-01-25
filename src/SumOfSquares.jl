module SumOfSquares

using LinearAlgebra

# MultivariatePolynomials extension

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets
export @set
using MultivariateMoments

include("matpoly.jl")
include("sosdec.jl")
include("certificate.jl")

# MOI extension

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

using PolyJuMP
export Poly

include("attributes.jl")
include("diagonally_dominant.jl")
include("sos_polynomial.jl")

# Bridges
const MOIB = MOI.Bridges
include("sos_polynomial_bridge.jl")
include("sos_polynomial_in_semialgebraic_set_bridge.jl")
include("diagonally_dominant_bridge.jl")
include("scaled_diagonally_dominant_bridge.jl")

# JuMP extension

import Reexport
Reexport.@reexport using JuMP

include("utilities.jl")
include("variable.jl")
include("constraint.jl")

function _add_bridges(model::JuMP.AbstractModel)
    JuMP.add_bridge(model, PolyJuMP.ZeroPolynomialBridge)
    JuMP.add_bridge(model, PolyJuMP.ZeroPolynomialInAlgebraicSetBridge)
    JuMP.add_bridge(model, PolyJuMP.PlusMinusBridge)
    JuMP.add_bridge(model, SOSPolynomialBridge)
    JuMP.add_bridge(model, SOSPolynomialInSemialgebraicSetBridge)
    JuMP.add_bridge(model, DiagonallyDominantBridge)
    JuMP.add_bridge(model, ScaledDiagonallyDominantBridge)
end

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SOSCone)
    PolyJuMP.setdefault!(data, PolyJuMP.PosDefPolyMatrix, SOSMatrixCone)
end

export SOSModel
function SOSModel(args...; kwargs...)
    model = Model(args...; kwargs...)
    _add_bridges(model)
    setpolymodule!(model, SumOfSquares)
    model
end

end # module
