module Bridges

import MathOptInterface as MOI
import SumOfSquares as SOS

include("Variable/Variable.jl")
include("Constraint/Constraint.jl")

function add_all_bridges(model, ::Type{T}) where {T}
    Variable.add_all_bridges(model, T)
    Constraint.add_all_bridges(model, T)
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{
        SOS.GramMatrixAttribute,
        SOS.MomentMatrixAttribute,
        SOS.SOSDecompositionAttribute,
    },
    bridge::MOI.Bridges.Constraint.VectorSlackBridge,
)
    return MOI.get(model, attr, bridge.slack_in_set)
end

# TODO bridges should redirect to `MOI.get_fallback` as well so that
# we can just use `Union{MOI.ConstraintIndex,MOI.Bridges.AbstractBridge}` in the `get_fallback` in `attributes.jl`
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.SOSDecompositionAttribute,
    bridge::Union{
        Variable.KernelBridge,
        Constraint.ImageBridge,
        Constraint.SOSPolynomialInSemialgebraicSetBridge,
    },
)
    gram = MOI.get(
        model,
        SOS.GramMatrixAttribute(;
            multiplier_index = attr.multiplier_index,
            result_index = attr.result_index,
        ),
        bridge,
    )
    return SOS.SOSDecomposition(gram, attr.ranktol, attr.dec)
end

end
