module Constraint

# Symmetric PSD matrix bridges
include("psd2x2_bridge.jl")
include("diagonally_dominant_bridge.jl")
include("scaled_diagonally_dominant_bridge.jl")

# SOS polynomial bridges
include("sos_polynomial_bridge.jl")
include("sos_polynomial_in_semialgebraic_set_bridge.jl")

end
