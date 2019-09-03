struct PositiveSemidefinite2x2Bridge{T} <: MOIB.Variable.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    rsoc::MOI.ConstraintIndex{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone}
end

function MOIB.Variable.bridge_constrained_variable(
    ::Type{PositiveSemidefinite2x2Bridge{T}}, model::MOI.ModelLike,
    s::SOS.PositiveSemidefinite2x2ConeTriangle) where {T}
    x, rsoc = MOI.add_constrained_variables(model, MOI.RotatedSecondOrderCone(3))
    return PositiveSemidefinite2x2Bridge{T}(x, rsoc)
end

function MOIB.Variable.supports_constrained_variable(
    ::Type{<:PositiveSemidefinite2x2Bridge}, ::Type{SOS.PositiveSemidefinite2x2ConeTriangle})
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:PositiveSemidefinite2x2Bridge})
    return [(MOI.RotatedSecondOrderCone,)]
end
function MOIB.added_constraint_types(::Type{PositiveSemidefinite2x2Bridge{T}}) where T
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

function _variable_map(i::MOIB.Variable.IndexInVector)
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
                   i::MOIB.Variable.IndexInVector)
    return bridge.variables[_variable_map(i)]
end

function MOI.get(model::MOI.ModelLike, attr::MOI.VariablePrimal,
                 bridge::PositiveSemidefinite2x2Bridge, i::MOIB.Variable.IndexInVector)
    value = MOI.get(model, attr, _variable(bridge, i))
    if i.value == 2
        value /= √2
    end
    return value
end

function MOIB.bridged_function(bridge::PositiveSemidefinite2x2Bridge{T},
                               i::MOIB.Variable.IndexInVector) where T
    func = MOI.SingleVariable(_variable(bridge, i))
    if i.value == 2
        return MOIU.operate(/, T, func, convert(T, √2))
    else
        return convert(MOI.ScalarAffineFunction{T}, func)
    end
end
function MOIB.Variable.unbridged_map(
    bridge::PositiveSemidefinite2x2Bridge{T},
    vi::MOI.VariableIndex, i::MOIB.Variable.IndexInVector) where T

    sv = MOI.SingleVariable(vi)
    if i.value == 2
        func = MOIU.operate(*, T, convert(T, √2), sv)
    else
        func = convert(MOI.ScalarAffineFunction{T}, sv)
    end
    return (_variable(bridge, i) => func,)
end
