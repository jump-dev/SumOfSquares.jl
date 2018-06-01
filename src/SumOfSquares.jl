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
import JuMP: validmodel, addtoexpr_reorder

include("variable.jl")
include("constraint.jl")

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.Poly{true}, SOSPoly)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, SOSCone)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPolyMatrix, SOSMatrixCone)
end

function SOSModel(; kwargs...)
    m = Model(; kwargs...)
    setpolymodule!(m, SumOfSquares)
    m
end

function PolyJuMP.getslack(c::SOSConstraint)
    getvalue(c.slack)
end

end # module
