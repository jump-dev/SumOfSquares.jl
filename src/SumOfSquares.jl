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

function SOSModel(args...; kwargs...)
    model = Model(args...; kwargs...)
    setpolymodule!(model, SumOfSquares)
    model
end

end # module
