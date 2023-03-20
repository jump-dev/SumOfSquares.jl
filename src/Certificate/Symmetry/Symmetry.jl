module Symmetry

import LinearAlgebra

import MutableArithmetics as MA
import MultivariatePolynomials as MP
import MultivariateBases as MB

import SymbolicWedderburn

import SumOfSquares

struct Pattern{GT,AT<:SymbolicWedderburn.Action}
    group::GT
    action::AT
end

include("utils.jl")
include("wedderburn.jl")
include("block_diag.jl")

end
