struct EmptyBridge{T} <: MOIB.AbstractBridge
end

function EmptyBridge{T}(model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
                        s::EmptyCone) where {T}
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    return EmptyBridge{T}()
end

function MOI.supports_constraint(::Type{<:EmptyBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:EmptyCone})
    return true
end
function MOIB.added_constraint_types(::Type{<:EmptyBridge{T}}) where {T}
    # TODO remove vov-in-Nonneg when MOI v0.8.3 is released
    return Tuple{DataType, DataType}[(MOI.VectorOfVariables, MOI.Nonnegatives)]
end
function MOIB.concrete_bridge_type(::Type{<:EmptyBridge{T}},
                                   ::Type{<:MOI.AbstractVectorFunction},
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
