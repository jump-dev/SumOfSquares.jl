struct EmptyBridge{T,F<:MOI.AbstractVectorFunction} <: MOI.Bridges.Constraint.AbstractBridge end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{EmptyBridge{T}},
    model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction,
    s::SOS.EmptyCone,
) where {T}
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    return EmptyBridge{T}()
end

function MOI.supports_constraint(
    ::Type{<:EmptyBridge},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.EmptyCone},
)
    return true
end
function MOI.Bridges.added_constrained_variable_types(::Type{<:EmptyBridge})
    return Tuple{DataType}[]
end
function MOI.Bridges.added_constraint_types(::Type{<:EmptyBridge})
    return Tuple{DataType,DataType}[]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:EmptyBridge{T}},
    ::Type{F},
    ::Type{SOS.EmptyCone},
) where {T,F<:MOI.AbstractVectorFunction}
    return EmptyBridge{T,F}
end

# Indices
function MOI.delete(::MOI.ModelLike, ::EmptyBridge) end

function MOI.get(
    ::MOI.ModelLike,
    ::ConstraintFunction,
    ::EmptyBridge{T,F},
) where {T,F}
    return MOI.Utilities.zero_with_output_dimension(F, 0)
end

function MOI.get(
    ::MOI.ModelLike,
    ::Union{MOI.ConstraintDual,MOI.ConstraintPrimal},
    bridge::EmptyBridge{T},
) where {T}
    return T[]
end
