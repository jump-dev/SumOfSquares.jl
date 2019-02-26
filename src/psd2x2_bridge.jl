# PSD constraints on 2x2 matrices are SOC representable.
# [Q11 Q12] is PSD iff Q11, Q22 ≥ 0 and       Q11*Q22 ≥     Q12 ^2
# [Q12 Q22]                             <=> 2*Q11*Q22 ≥ (√2*Q12)^2
struct PositiveSemidefinite2x2Bridge{T, F} <: MOIB.AbstractBridge
    rsoc::MOI.ConstraintIndex{F, MOI.RotatedSecondOrderCone}
end

function PositiveSemidefinite2x2Bridge{T, F}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::PositiveSemidefinite2x2ConeTriangle) where {T, F}
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    fs = MOIU.eachscalar(f)
    g = MOIU.operate(vcat, T, fs[1], fs[3], √2 * fs[2])
    rsoc = MOI.add_constraint(model, g, MOI.RotatedSecondOrderCone(3))
    return PositiveSemidefinite2x2Bridge{T, F}(rsoc)
end

function MOI.supports_constraint(::Type{<:PositiveSemidefinite2x2Bridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{PositiveSemidefinite2x2ConeTriangle})
    return true
end
function MOIB.added_constraint_types(::Type{PositiveSemidefinite2x2Bridge{T, F}}) where {T, F}
    return [(F, MOI.RotatedSecondOrderCone)]
end
function MOIB.concrete_bridge_type(
    ::Type{<:PositiveSemidefinite2x2Bridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{PositiveSemidefinite2x2ConeTriangle}) where T
    S = MOIU.scalar_type(F)
    G = MOIU.promote_operation(*, T, T, S)
    H = MOIU.promote_operation(vcat, T, G)
    return PositiveSemidefinite2x2Bridge{T, H}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::PositiveSemidefinite2x2Bridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.RotatedSecondOrderCone}) where {T, F}
    return 1
end
function MOI.get(bridge::PositiveSemidefinite2x2Bridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.RotatedSecondOrderCone}) where {T, F}
    return [bridge.rsoc]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::PositiveSemidefinite2x2Bridge)
    MOI.delete(model, bridge.rsoc)
end

# TODO ConstraintPrimal
function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::PositiveSemidefinite2x2Bridge)
    dual = MOI.get(model, attr, bridge.rsoc)
    return [dual[1], dual[3] / √2, dual[2]]
end
