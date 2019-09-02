struct SOS.PositiveSemidefinite2x2Bridge{T} <: MOIB.Variable.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    rsoc::MOI.ConstraintIndex{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone}
end

function MOIB.Variable.bridge_constrained_variable(
    ::Type{SOS.PositiveSemidefinite2x2VariableBridge{T}}, model::MOI.ModelLike,
    s::SOS.PositiveSemidefinite2x2ConeTriangle) where {T}
    x, rsoc = MOI.add_constrained_variables(model, MOI.RotatedSecondOrderCone(3))
    Q12 = MOIU.operate(/, T, MOI.SingleVariable(x[3]), convert(T, √2))
    g = typeof(Q12)[MOI.SingleVariable(x[1]), Q12, MOI.SingleVariable(x[2])]
    return SOS.PositiveSemidefinite2x2VariableBridge{T}(x, rsoc)
end

function MOIB.Variable.supports_constrained_variable(
    ::Type{<:SOS.PositiveSemidefinite2x2VariableBridge}, ::Type{SOS.PositiveSemidefinite2x2ConeTriangle})
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:SOS.PositiveSemidefinite2x2VariableBridge})
    return [(MOI.RotatedSecondOrderCone,)]
end
function MOIB.added_constraint_types(::Type{RSOCtoPSDBridge{T}}) where T
    return Tuple{DataType, DataType}[]
end

# Attributes, Bridge acting as a model
function MOI.get(::SOS.PositiveSemidefinite2x2VariableBridge,
                 ::MOI.NumberOfVariables)
    return 3
end
function MOI.get(bridge::SOS.PositiveSemidefinite2x2VariableBridge,
                 ::MOI.ListOfVariableIndices)
    return bridge.variables
end
function MOI.get(
    ::SOS.PositiveSemidefinite2x2VariableBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone})
    return 1
end
function MOI.get(
    bridge::SOS.PositiveSemidefinite2x2VariableBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,
                                  MOI.RotatedSecondOrderCone})
    return [bridge.rsoc]
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::SOS.PositiveSemidefinite2x2VariableBridge)
    MOI.delete(model, bridge.variables)
end

# Attributes, Bridge acting as a constraint

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet,
                 ::SOS.PositiveSemidefinite2x2VariableBridge)
    return SOS.PositiveSemidefinite2x2ConeTriangle()
end

function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintPrimal,
                 bridge::SOS.PositiveSemidefinite2x2VariableBridge)
    value = MOI.get(model, attr, bridge.rsoc)
    return [value[1], value[3] / √2, value[2]]
end
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::SOS.PositiveSemidefinite2x2VariableBridge)
    dual = MOI.get(model, attr, bridge.rsoc)
    # / 2 (because of different scalar product) * √2 (because of A^{-*} = 1 / √2
    return [dual[1], dual[3] / √2, dual[2]]
end

function _variable_map(i::IndexInVector)
    if i.value == 1
        return 1
    elseif i.value == 2
        return 3
    else
        @assert i.value == 3
        return 2
    end
end
function _variable(bridge::SOS.PositiveSemidefinite2x2VariableBridge,
                   i::IndexInVector)
    return bridge.variables[_variable_map(i)]
end

function MOI.get(model::MOI.ModelLike, attr::MOI.VariablePrimal,
                 bridge::SOS.PositiveSemidefinite2x2VariableBridge, i::IndexInVector)
    value = MOI.get(model, attr, _variable(bridge, i))
    if i.value == 2
        value /= √2
    end
    return value
end

function MOIB.bridged_function(bridge::SOS.PositiveSemidefinite2x2VariableBridge{T},
                               i::IndexInVector) where T
    func = MOI.SingleVariable(_variable(bridge, i))
    if i.value == 2
        return MOIU.operate(/, T, func, convert(T, √2))
    else
        return convert(MOI.ScalarAffineFunction{T}, func)
    end
end
function MOIB.Variable.unbridged_map(
    bridge::SOS.PositiveSemidefinite2x2VariableBridge{T},
    vi::MOI.VariableIndex, i::IndexInVector) where T

    sv = MOI.SingleVariable(vi)
    if i.value == 2
        func = MOIU.operate(*, T, convert(T, √2), sv)
    else
        func = convert(MOI.ScalarAffineFunction{T}, sv)
    end
    return (_variable(bridge, i) => func,)
end
