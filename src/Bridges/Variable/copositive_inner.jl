struct CopositiveInnerBridge{T,S} <: MOI.Bridges.Variable.AbstractBridge
    matrix_variables::Vector{MOI.VariableIndex}
    matrix_constraint::MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    nonneg_variables::Vector{MOI.VariableIndex}
    nonneg_constraint::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        MOI.Nonnegatives,
    }
end

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{CopositiveInnerBridge{T,S}},
    model::MOI.ModelLike,
    set::SOS.CopositiveInnerCone,
) where {T,S}
    side_dimension = MOI.side_dimension(set.psd_inner)
    num_off_diag =
        MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(side_dimension - 1))
    return CopositiveInnerBridge{T,S}(
        MOI.add_constrained_variables(model, set.psd_inner)...,
        MOI.add_constrained_variables(model, MOI.Nonnegatives(num_off_diag))...,
    )
end

function MOI.Bridges.Variable.supports_constrained_variable(
    ::Type{<:CopositiveInnerBridge},
    ::Type{<:SOS.CopositiveInnerCone},
)
    return true
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{CopositiveInnerBridge{T,S}},
) where {T,S}
    return [(S,), (MOI.Nonnegatives,)]
end
function MOI.Bridges.added_constraint_types(::Type{<:CopositiveInnerBridge})
    return Tuple{Type,Type}[]
end
function MOI.Bridges.Variable.concrete_bridge_type(
    ::Type{<:CopositiveInnerBridge{T}},
    ::Type{SOS.CopositiveInnerCone{S}},
) where {T,S}
    return CopositiveInnerBridge{T,S}
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::CopositiveInnerBridge, ::MOI.NumberOfVariables)
    return length(bridge.matrix_variables) + length(bridge.nonneg_variables)
end
function MOI.get(bridge::CopositiveInnerBridge, ::MOI.ListOfVariableIndices)
    return Iterators.flatten((bridge.matrix_variables, bridge.nonneg_variables))
end
function MOI.get(
    bridge::CopositiveInnerBridge{T,S},
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
) where {T,S}
    return 1
end
function MOI.get(
    bridge::CopositiveInnerBridge{T,S},
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {T,S}
    return [bridge.matrix_constraint]
end
function MOI.get(
    bridge::CopositiveInnerBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,MOI.Nonnegatives},
)
    return 1
end
function MOI.get(
    bridge::CopositiveInnerBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.Nonnegatives},
)
    return [bridge.nonneg_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::CopositiveInnerBridge)
    MOI.delete(model, bridge.matrix_variables)
    return MOI.delete(model, bridge.nonneg_variables)
end

# Attributes, Bridge acting as a constraint

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintSet,
    bridge::CopositiveInnerBridge,
)
    return SOS.CopositiveInnerCone(
        MOI.get(model, attr, bridge.matrix_constraint),
    )
end

# TODO ConstraintPrimal, ConstraintDual

# Vector index for the vectorization of the off-diagonal triangular part.
function offdiag_vector_index(i, j)
    if i < j
        return MOI.Utilities.trimap(i, j - 1)
    else
        throw(ArgumentError("Not off-diagonal"))
    end
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::CopositiveInnerBridge,
    i::MOI.Bridges.IndexInVector,
)
    value = MOI.get(model, attr, bridge.matrix_variables[i.value])
    row, col = MOI.Utilities.inverse_trimap(i.value)
    if row != col
        value += MOI.get(
            model,
            attr,
            bridge.nonneg_variables[offdiag_vector_index(i, j)],
        )
    end
    return value
end

function MOI.Bridges.bridged_function(
    bridge::CopositiveInnerBridge{T},
    i::MOI.Bridges.IndexInVector,
) where {T}
    func =
        convert(MOI.ScalarAffineFunction{T}, bridge.matrix_variables[i.value])
    row, col = MOI.Utilities.inverse_trimap(i.value)
    if row != col
        func = MOI.Utilities.operate!(
            +,
            T,
            func,
            bridge.nonneg_variables[MOI.Utilities.trimap(row, col - 1)],
        )
    end
    return func
end
function MOI.Bridges.Variable.unbridged_map(
    bridge::CopositiveInnerBridge{T},
    vi::MOI.VariableIndex,
    i::MOI.Bridges.IndexInVector,
) where {T}
    F = MOI.ScalarAffineFunction{T}
    func = convert(F, vi)
    map = bridge.matrix_variables[i.value] => func
    row, col = MOI.Utilities.inverse_trimap(i.value)
    if row == col
        return (map,)
    else
        nneg = bridge.nonneg_variables[MOI.Utilities.trimap(row, col - 1)]
        return (map, nneg => zero(F))
    end
end
