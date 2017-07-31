__precompile__()

module SumOfSquares

export SOSModel

import Base.show, Base.length, Base.getindex, Base.vect, Base.isless, Base.isempty, Base.start, Base.done, Base.next, Base.convert, Base.dot

using MultivariatePolynomials
using Polyhedra

include("certificate.jl")

using PolyJuMP, JuMP
import JuMP: validmodel, addtoexpr_reorder

include("variable.jl")
include("constraint.jl")

function SOSModel(; kwargs...)
    m = Model(; kwargs...)
    setpolymodule!(m, SumOfSquares)
    m
end

function PolyJuMP.getslack(c::SOSConstraint)
    getvalue(c.slack)
end

end # module
