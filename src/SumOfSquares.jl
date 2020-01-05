module SumOfSquares

using LinearAlgebra

import Reexport

import MutableArithmetics
const MA = MutableArithmetics

# MultivariatePolynomials extension

import MultivariatePolynomials
const MP = MultivariatePolynomials
# @set assumes that `SemialgebraicSets` is defined
Reexport.@reexport using SemialgebraicSets
Reexport.@reexport using MultivariateMoments
_promote_sum(T::Type) = MA.promote_operation(+, T, T)
_promote_add_mul(T::Type) = MA.promote_operation(MA.add_mul, T, T, T)

include("gram_matrix.jl")

using PolyJuMP
abstract type SOSLikeCone <: PolyJuMP.PolynomialSet end
Base.broadcastable(cone::SOSLikeCone) = Ref(cone)
# FIXME Maybe replace by PSDLike MOI sets

function matrix_cone_type end

include("Certificate/Certificate.jl")
include("rand.jl")

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
include("Bridges/Bridges.jl")

# JuMP extension

Reexport.@reexport using JuMP

include("sosdec.jl")
include("utilities.jl")
include("variable.jl")
include("constraint.jl")

# Graphs 
include("csp/ChordalExtensionGraph.jl")
const CEG = SumOfSquares.ChordalExtensionGraph
const APL = MP.AbstractPolynomialLike
include("csp/csp_graph.jl")
include("csp/sparse_putinar.jl")


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
