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

using PolyJuMP
export Poly

include("sos_polynomial.jl")

# Bridges
const MOIB = MOI.Bridges
include("sos_polynomial_bridge.jl")
include("sos_polynomial_in_semialgebraic_set_bridge.jl")

# JuMP extension

import Reexport
Reexport.@reexport using JuMP

include("variable.jl")
include("constraint.jl")

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SOSCone)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPolyMatrix, SOSMatrixCone)
end

export SOSModel
function SOSModel(args...; kwargs...)
    model = Model(args...; kwargs...)
    setpolymodule!(model, SumOfSquares)
    model
end

end # module
