module Symmetry

import LinearAlgebra

import MutableArithmetics
const MA = MutableArithmetics
import MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases
const MB = MultivariateBases

import SymbolicWedderburn

import SumOfSquares

struct Pattern{GT,AT<:SymbolicWedderburn.Action}
    group::GT
    action::AT
end

include("wedderburn.jl")
include("block_diag.jl")

end
