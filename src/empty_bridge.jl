struct EmptyBridge{T} <: MOIB.Constraint.AbstractBridge
end

function MOIB.Constraint.bridge_constraint(
    ::Type{EmptyBridge{T}}, model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::EmptyCone) where {T}
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    return EmptyBridge{T}()
end

function MOI.supports_constraint(::Type{<:EmptyBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:EmptyCone})
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:EmptyBridge})
    return Tuple{DataType}[]
end
function MOIB.added_constraint_types(::Type{<:EmptyBridge})
    return Tuple{DataType, DataType}[]
end
function MOIB.Constraint.concrete_bridge_type(
    ::Type{<:EmptyBridge{T}}, ::Type{<:MOI.AbstractVectorFunction},
    ::Type{EmptyCone}) where T
    return EmptyBridge{T}
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::EmptyBridge) end

# TODO ConstraintPrimal
function MOI.get(::MOI.ModelLike,
                 ::Union{MOI.ConstraintDual, MOI.ConstraintPrimal},
                 bridge::EmptyBridge{T}) where T
    return T[]
end
