module Constraint

# Symmetric PSD matrix bridges
include("empty.jl")
include("psd2x2.jl")
include("diagonally_dominant.jl")
include("scaled_diagonally_dominant.jl")

# SOS polynomial bridges
include("sos_polynomial.jl")
include("sos_polynomial_in_semialgebraic_set.jl")

end
