struct EmptyBridge{T,F<:MOI.AbstractVectorFunction} <:
       MOI.Bridges.Constraint.AbstractBridge end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{EmptyBridge{T,F}},
    model::MOI.ModelLike,
    f::F,
    s::SOS.EmptyCone,
) where {T,F}
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    return EmptyBridge{T,F}()
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

function _empty_function(::Type{MOI.VectorOfVariables})
    return MOI.VectorOfVariables(MOI.VariableIndex[])
end
function _empty_function(::Type{F}) where {F}
    return MOI.Utilities.zero_with_output_dimension(F, 0)
end
function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    ::EmptyBridge{T,F},
) where {T,F}
    return _empty_function(F)
end

function MOI.get(
    ::MOI.ModelLike,
    ::Union{MOI.ConstraintDual,MOI.ConstraintPrimal},
    bridge::EmptyBridge{T},
) where {T}
    return T[]
end
