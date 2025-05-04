module Variable

import SparseArrays
import MutableArithmetics as MA
import StarAlgebras as SA
import MultivariateBases as MB
import MathOptInterface as MOI
import MultivariatePolynomials as MP
import SumOfSquares as SOS

include("psd2x2.jl")
include("scaled_diagonally_dominant.jl")
include("copositive_inner.jl")
include("kernel.jl")

function add_all_bridges(model, ::Type{T}) where {T}
    MOI.Bridges.add_bridge(model, PositiveSemidefinite2x2Bridge{T})
    MOI.Bridges.add_bridge(model, ScaledDiagonallyDominantBridge{T})
    MOI.Bridges.add_bridge(model, CopositiveInnerBridge{T})
    MOI.Bridges.add_bridge(model, KernelBridge{T})
    MOI.Bridges.add_bridge(model, LowRankBridge{T})
    return
end

end
