module SumOfSquares

using LinearAlgebra

import Reexport

# MultivariatePolynomials extension

import MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets
export @set
Reexport.@reexport using MultivariateMoments

include("gram_matrix.jl")
include("sosdec.jl")
include("certificate.jl")

# MOI extension

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

Reexport.@reexport using PolyJuMP

include("attributes.jl")
include("psd2x2.jl")
include("diagonally_dominant.jl")
include("sos_polynomial.jl")
include("copositive_inner.jl")

# Bridges
const MOIB = MOI.Bridges

# Variable Bridges
include("variable_bridge.jl")
include("scaled_diagonally_dominant_variable_bridge.jl")
include("generic_variable_bridge.jl")
include("copositive_inner_variable_bridge.jl")

# Constraint Bridges
include("sos_polynomial_bridge.jl")
include("sos_polynomial_in_semialgebraic_set_bridge.jl")
include("diagonally_dominant_bridge.jl")
include("empty_bridge.jl")
include("psd2x2_bridge.jl")
include("scaled_diagonally_dominant_bridge.jl")

# JuMP extension

Reexport.@reexport using JuMP

include("utilities.jl")
include("variable.jl")
include("constraint.jl")

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SOSCone)
    PolyJuMP.setdefault!(data, PolyJuMP.PosDefPolyMatrix, SOSMatrixCone)
    return
end

export SOSModel
function SOSModel(args...; kwargs...)
    model = Model(args...; kwargs...)
    setpolymodule!(model, SumOfSquares)
    model
end

end # module
