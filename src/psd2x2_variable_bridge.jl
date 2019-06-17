struct PositiveSemidefinite2x2VariableBridge{T} <: AbstractVariableBridge
    variables::Vector{MOI.VariableIndex}
    rsoc::MOI.ConstraintIndex{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone}
end

function add_variable_bridge(
    ::Type{PositiveSemidefinite2x2VariableBridge{T}}, model::MOI.ModelLike,
    s::PositiveSemidefinite2x2ConeTriangle) where {T}
    x = MOI.add_variables(model, 3)
    rsoc = MOI.add_constraint(model, MOI.VectorOfVariables(x), MOI.RotatedSecondOrderCone(3))
    Q12 = MOIU.operate(/, T, MOI.SingleVariable(x[3]), convert(T, √2))
    g = typeof(Q12)[MOI.SingleVariable(x[1]), Q12, MOI.SingleVariable(x[2])]
    return g, PositiveSemidefinite2x2VariableBridge{T}(x, rsoc)
end

function MOIB.added_constraint_types(
    ::Type{PositiveSemidefinite2x2VariableBridge{T}}) where {T}
    return [(MOI.VectorOfVariables, MOI.RotatedSecondOrderCone)]
end

function variable_bridge_type(::Type{PositiveSemidefinite2x2ConeTriangle},
                              T::Type)
    return PositiveSemidefinite2x2VariableBridge{T}
end


# Attributes, VariableBridge acting as an model
function MOI.get(::PositiveSemidefinite2x2VariableBridge,
                 ::MOI.NumberOfVariables)
    return 3
end
function MOI.get(
    ::PositiveSemidefinite2x2VariableBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,
                              MOI.RotatedSecondOrderCone})
    return 1
end
function MOI.get(
    bridge::PositiveSemidefinite2x2VariableBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,
                                  MOI.RotatedSecondOrderCone})
    return [bridge.rsoc]
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::PositiveSemidefinite2x2VariableBridge)
    MOI.delete(model, bridge.rsoc)
    for vi in bridge.variables
        MOI.delete(model, vi)
    end
end

function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintPrimal,
                 bridge::PositiveSemidefinite2x2VariableBridge)
    value = MOI.get(model, attr, bridge.rsoc)
    return [value[1], value[3] / √2, value[2]]
end
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::PositiveSemidefinite2x2VariableBridge)
    dual = MOI.get(model, attr, bridge.rsoc)
    # / 2 (because of different scalar product) * √2 (because of A^{-*} = 1 / √2
    return [dual[1], dual[3] / √2, dual[2]]
end
