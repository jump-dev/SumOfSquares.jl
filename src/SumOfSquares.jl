module SumOfSquares

using LinearAlgebra

import Reexport

import MutableArithmetics
const MA = MutableArithmetics

# MultivariatePolynomials extension

import MultivariatePolynomials
const MP = MultivariatePolynomials
Reexport.@reexport using MultivariateBases
# @set assumes that `SemialgebraicSets` is defined
Reexport.@reexport using SemialgebraicSets
Reexport.@reexport using MultivariateMoments

include("gram_matrix.jl")
include("sosdec.jl")

using PolyJuMP
abstract type SOSLikeCone <: PolyJuMP.PolynomialSet end
Base.broadcastable(cone::SOSLikeCone) = Ref(cone)
# FIXME Maybe replace by PSDLike MOI sets

function matrix_cone_type end

export Certificate
include("Certificate/Certificate.jl")
using .Certificate: Sparsity, SignSymmetry, ChordalCompletion, ClusterCompletion, Symmetry
export Sparsity, SignSymmetry, ChordalCompletion, ClusterCompletion
export Symmetry

# TODO remove in SumOfSquares v1.0
@deprecate NoSparsity Sparsity.NoPattern
@deprecate VariableSparsity Sparsity.Variable
@deprecate MonomialSparsity Sparsity.Monomial

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
