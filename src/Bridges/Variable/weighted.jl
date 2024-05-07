struct KernelBridge{T,M} <: MOI.Bridges.Variable.AbstractBridge
    variables::Vector{Vector{MOI.VariableIndex}}
    constraints::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables}}
end

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{KernelBridge{T}},
    model::MOI.ModelLike,
    set::SOS.WeightedSOSCone,
) where {T}
    variables = Vector{MOI.VariableIndex}[]
    constraints = MOI.ConstraintIndex{MOI.VectorOfVariables}[]
    acc = MA.Zero()
    for (gram_basis, weight) in zip(set.gram_bases, set.weights)
        gram, vars, con = SOS.add_gram_matrix(model, M, gram_basis, T)
        push!(variables, vars)
        push!(constraints, con)
        acc = MA.add_mul!!(acc, weight, gram)
    end
    return KernelBridge{T}(variables, constraints)
end

function MOI.Bridges.Variable.supports_constrained_variable(
    ::Type{<:KernelBridge},
    ::Type{<:SOS.WeightedSOSCone},
)
    return true
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{KernelBridge{T,M}},
) where {T,M}
    return constrained_variable_types(M)
end
function MOI.Bridges.added_constraint_types(
    ::Type{PositiveSemidefinite2x2Bridge{T}},
) where {T}
    return Tuple{Type,Type}[]
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::KernelBridge, ::MOI.NumberOfVariables)
    return sum(length, bridge.variables)
end
function MOI.get(
    bridge::KernelBridge,
    ::MOI.ListOfVariableIndices,
)
    return reduce(vcat, bridge.variables)
end
function MOI.get(
    ::KernelBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
) where {S<:MOI.AbstractVectorSet}
    return count(bridge.constraints) do ci
        ci isa MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    end
end
function MOI.get(
    bridge::KernelBridge,
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        S,
    },
)
    return [ci for ci in bridge.constraints if ci isa MOI.ConstraintIndex{MOI.VectorOfVariables,S}]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::KernelBridge)
    for vars in bridge.variables
        MOI.delete(model, vars)
    end
    return
end

# Attributes, Bridge acting as a constraint

#function MOI.get(
#    ::MOI.ModelLike,
#    ::MOI.ConstraintSet,
#    ::KernelBridge,
#)
#    return SOS.WeightedSOSCone()
#end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintPrimal,
    bridge::KernelBridge,
)
    value = MOI.get(model, attr, bridge.rsoc)
    return [value[1], value[3] / √2, value[2]]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::KernelBridge,
)
    dual = MOI.get(model, attr, bridge.rsoc)
    # / 2 (because of different scalar product) * √2 (because of A^{-*} = 1 / √2
    return [dual[1], dual[3] / √2, dual[2]]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::KernelBridge,
    i::MOI.Bridges.IndexInVector,
)
    value = MOI.get(model, attr, _variable(bridge, i))
    if i.value == 2
        value /= √2
    end
    return value
end

function MOI.Bridges.bridged_function(
    bridge::KernelBridge{T},
    i::MOI.Bridges.IndexInVector,
) where {T}
    func = _variable(bridge, i)
    if i.value == 2
        return MOI.Utilities.operate(/, T, func, convert(T, √2))
    else
        return convert(MOI.ScalarAffineFunction{T}, func)
    end
end
