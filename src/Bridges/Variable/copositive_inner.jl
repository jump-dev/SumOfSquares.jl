struct CopositiveInnerBridge{T, S} <: MOIB.Variable.AbstractBridge
    matrix_variables::Vector{MOI.VariableIndex}
    matrix_constraint::MOI.ConstraintIndex{MOI.VectorOfVariables, S}
    nonneg_variables::Vector{MOI.VariableIndex}
    nonneg_constraint::MOI.ConstraintIndex{MOI.VectorOfVariables, MOI.Nonnegatives}
end

function MOIB.Variable.bridge_constrained_variable(
    ::Type{CopositiveInnerBridge{T, S}},
    model::MOI.ModelLike, set::SOS.CopositiveInnerCone) where {T, S}

    side_dimension = MOI.side_dimension(set.psd_inner)
    num_off_diag = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(side_dimension - 1))
    return CopositiveInnerBridge{T, S}(
        MOI.add_constrained_variables(model, set.psd_inner)...,
        MOI.add_constrained_variables(model, MOI.Nonnegatives(num_off_diag))...
    )
end

function MOIB.Variable.supports_constrained_variable(
    ::Type{<:CopositiveInnerBridge}, ::Type{<:SOS.CopositiveInnerCone})
    return true
end
function MOIB.added_constrained_variable_types(::Type{CopositiveInnerBridge{T, S}}) where {T, S}
    return [(S,), (MOI.Nonnegatives,)]
end
function MOIB.added_constraint_types(::Type{<:CopositiveInnerBridge})
    return Tuple{DataType, DataType}[]
end
function MOIB.Variable.concrete_bridge_type(
    ::Type{<:CopositiveInnerBridge{T}}, ::Type{SOS.CopositiveInnerCone{S}}) where {T, S}
    return CopositiveInnerBridge{T, S}
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::CopositiveInnerBridge,
                 ::MOI.NumberOfVariables)
    return length(bridge.matrix_variables) + length(bridge.nonneg_variables)
end
function MOI.get(bridge::CopositiveInnerBridge,
                 ::MOI.ListOfVariableIndices)
    return Iterators.flatten((bridge.matrix_variables, bridge.nonneg_variables))
end
function MOI.get(bridge::CopositiveInnerBridge{T, S},
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, S}) where {T, S}
    return 1
end
function MOI.get(bridge::CopositiveInnerBridge{T, S},
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}) where {T, S}
    return [bridge.matrix_constraint]
end
function MOI.get(bridge::CopositiveInnerBridge,
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, MOI.Nonnegatives})
    return 1
end
function MOI.get(bridge::CopositiveInnerBridge,
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, MOI.Nonnegatives})
    return [bridge.nonneg_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::CopositiveInnerBridge)
    MOI.delete(model, bridge.matrix_variables)
    MOI.delete(model, bridge.nonneg_variables)
end

# Attributes, Bridge acting as a constraint

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintSet,
                 bridge::CopositiveInnerBridge)
    return SOS.CopositiveInnerCone(MOI.get(model, attr, bridge.matrix_constraint))
end

# TODO ConstraintPrimal, ConstraintDual

# See https://www.juliaopt.org/MathOptInterface.jl/v0.9.1/apireference/#MathOptInterface.AbstractSymmetricMatrixSetTriangle
function matrix_indices(k)
    j = div(1 + isqrt(8k - 7), 2)
    i = k - div((j - 1) * j, 2)
    return i, j
end
# Vector index for the vectorization of the triangular part.
function vector_index(i, j)
    return div((j - 1) * j, 2) + i
end
# Vector index for the vectorization of the off-diagonal triangular part.
function offdiag_vector_index(i, j)
    if i < j
        return vector_index(i, j - 1)
    else
        throw(ArgumentError())
    end
end

function MOI.get(model::MOI.ModelLike, attr::MOI.VariablePrimal,
                 bridge::CopositiveInnerBridge, i::MOIB.Variable.IndexInVector)
    value = MOI.get(model, attr, bridge.matrix_variables[i.value])
    row, col = matrix_indices(i.value)
    if row != col
        value += MOI.get(model, attr, bridge.nonneg_variables[offdiag_vector_index(i, j)])
    end
    return value
end

function MOIB.bridged_function(bridge::CopositiveInnerBridge{T},
                               i::MOIB.Variable.IndexInVector) where T
    func = convert(MOI.ScalarAffineFunction{T},
                   MOI.SingleVariable(bridge.matrix_variables[i.value]))
    row, col = matrix_indices(i.value)
    if row != col
        func = MOIU.operate!(+, T, func, MOI.SingleVariable(
            bridge.nonneg_variables[vector_index(row, col - 1)]))
    end
    return func
end
function MOIB.Variable.unbridged_map(
    bridge::CopositiveInnerBridge{T},
    vi::MOI.VariableIndex, i::MOIB.Variable.IndexInVector) where T

    # TODO
    return nothing
end
