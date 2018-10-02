__precompile__()

module SumOfSquares

export SOSModel

using Compat, Compat.LinearAlgebra

using MultivariatePolynomials
const MP = MultivariatePolynomials
using SemialgebraicSets
export @set
using MultivariateMoments

include("matpoly.jl")
include("sosdec.jl")

include("certificate.jl")

using JuMP
using PolyJuMP
export Poly

include("variable.jl")
include("constraint.jl")

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SOSCone)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPolyMatrix, SOSMatrixCone)
end

function SOSModel(; kwargs...)
    m = Model(; kwargs...)
    setpolymodule!(m, SumOfSquares)
    m
end

end # module
