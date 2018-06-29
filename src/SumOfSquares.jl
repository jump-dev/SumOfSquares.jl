__precompile__()

module SumOfSquares

export SOSModel

using MultivariatePolynomials
const MP = MultivariatePolynomials
using MultivariateMoments
using SemialgebraicSets

include("matpoly.jl")
include("sosdec.jl")

include("certificate.jl")

using PolyJuMP, JuMP

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
