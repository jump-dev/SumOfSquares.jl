struct GenericVariableBridge{T, S} <: AbstractVariableBridge
    gram_matrix::Vector{MOI.VariableIndex}
    gram_constraint::MOI.ConstraintIndex{MOI.VectorOfVariables, S}
end

function add_variable_bridge(::Type{GenericVariableBridge{T, S}},
                             model::MOI.ModelLike, set::S) where {T, S}
    gram_matrix = MOI.add_variables(model, MOI.dimension(set))
    ci = MOI.add_constraint(model, MOI.VectorOfVariables(gram_matrix), set)
    func = MOI.SingleVariable[MOI.SingleVariable(v) for v in gram_matrix]
    # Need to specify `T` in the constructor as it cannot be deduced
    return func, GenericVariableBridge{T, S}(gram_matrix, ci)
end

function MOIB.added_constraint_types(::Type{GenericVariableBridge{T, S}}) where {T, S}
    return [(MOI.VectorOfVariables, S)]
end

function variable_bridge_type(S::Type{<:MOI.AbstractVectorSet}, T::Type)
    return GenericVariableBridge{T, S}
end

# Attributes, VariableBridge acting as an model
function MOI.get(bridge::GenericVariableBridge, ::MOI.NumberOfVariables)
    return length(bridge.gram_matrix)
end
function MOI.get(bridge::GenericVariableBridge{T, S},
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, S}) where {T, S}
    return bridge.gram_constraint isa MOI.ConstraintIndex{MOI.VectorOfVariables, S} ? 1 : 0
end
function MOI.get(bridge::GenericVariableBridge{T, S},
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}) where {T, S}
    return [bridge.gram_constraint]
end

function MOI.delete(model::MOI.ModelLike, bridge::GenericVariableBridge)
    # First delete the constraint in which the Gram matrix appears
    MOI.delete(model, bridge.gram_constraint)
    # Now we delete the Gram matrix
    for variable in bridge.gram_matrix
        MOI.delete(model, variable)
    end
end

function MOI.get(model::MOI.ModelLike, ::MomentMatrixAttribute,
                 bridge::GenericVariableBridge)
    return MOI.get(model, MOI.ConstraintDual(), bridge.gram_constraint)
end

function MOI.get(model::MOI.ModelLike, ::GramMatrixAttribute,
                 bridge::GenericVariableBridge)
    return MOI.get(model, MOI.VariablePrimal(), bridge.gram_matrix)
end
