struct PositiveSemidefinite2x2Bridge{T} <: MOI.Bridges.Variable.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    rsoc::MOI.ConstraintIndex{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone}
end

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{PositiveSemidefinite2x2Bridge{T}}, model::MOI.ModelLike,
    s::SOS.PositiveSemidefinite2x2ConeTriangle) where {T}
    x, rsoc = MOI.add_constrained_variables(model, MOI.RotatedSecondOrderCone(3))
    return PositiveSemidefinite2x2Bridge{T}(x, rsoc)
end

function MOI.Bridges.Variable.supports_constrained_variable(
    ::Type{<:PositiveSemidefinite2x2Bridge}, ::Type{SOS.PositiveSemidefinite2x2ConeTriangle})
    return true
end
function MOI.Bridges.added_constrained_variable_types(::Type{<:PositiveSemidefinite2x2Bridge})
    return [(MOI.RotatedSecondOrderCone,)]
end
function MOI.Bridges.added_constraint_types(::Type{PositiveSemidefinite2x2Bridge{T}}) where T
    return Tuple{DataType, DataType}[]
end

# Attributes, Bridge acting as a model
function MOI.get(::PositiveSemidefinite2x2Bridge,
                 ::MOI.NumberOfVariables)
    return 3
end
function MOI.get(bridge::PositiveSemidefinite2x2Bridge,
                 ::MOI.ListOfVariableIndices)
    return bridge.variables
end
function MOI.get(
    ::PositiveSemidefinite2x2Bridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone})
    return 1
end
function MOI.get(
    bridge::PositiveSemidefinite2x2Bridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,
                                  MOI.RotatedSecondOrderCone})
    return [bridge.rsoc]
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::PositiveSemidefinite2x2Bridge)
    MOI.delete(model, bridge.variables)
end

# Attributes, Bridge acting as a constraint

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet,
                 ::PositiveSemidefinite2x2Bridge)
    return SOS.PositiveSemidefinite2x2ConeTriangle()
end

function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintPrimal,
                 bridge::PositiveSemidefinite2x2Bridge)
    value = MOI.get(model, attr, bridge.rsoc)
    return [value[1], value[3] / √2, value[2]]
end
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::PositiveSemidefinite2x2Bridge)
    dual = MOI.get(model, attr, bridge.rsoc)
    # / 2 (because of different scalar product) * √2 (because of A^{-*} = 1 / √2
    return [dual[1], dual[3] / √2, dual[2]]
end

function _variable_map(i::MOI.Bridges.IndexInVector)
    if i.value == 1
        return 1
    elseif i.value == 2
        return 3
    else
        @assert i.value == 3
        return 2
    end
end
function _variable(bridge::PositiveSemidefinite2x2Bridge,
                   i::MOI.Bridges.IndexInVector)
    return bridge.variables[_variable_map(i)]
end

function MOI.get(model::MOI.ModelLike, attr::MOI.VariablePrimal,
                 bridge::PositiveSemidefinite2x2Bridge, i::MOI.Bridges.IndexInVector)
    value = MOI.get(model, attr, _variable(bridge, i))
    if i.value == 2
        value /= √2
    end
    return value
end

function MOI.Bridges.bridged_function(bridge::PositiveSemidefinite2x2Bridge{T},
                               i::MOI.Bridges.IndexInVector) where T
    func = _variable(bridge, i)
    if i.value == 2
        return MOI.Utilities.operate(/, T, func, convert(T, √2))
    else
        return convert(MOI.ScalarAffineFunction{T}, func)
    end
end
function MOI.Bridges.Variable.unbridged_map(
    bridge::PositiveSemidefinite2x2Bridge{T},
    vi::MOI.VariableIndex, i::MOI.Bridges.IndexInVector) where T

    if i.value == 2
        func = MOI.Utilities.operate(*, T, convert(T, √2), vi)
    else
        func = convert(MOI.ScalarAffineFunction{T}, vi)
    end
    return (_variable(bridge, i) => func,)
end
