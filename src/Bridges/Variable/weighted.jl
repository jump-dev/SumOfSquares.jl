struct KernelBridge{T,M} <: MOI.Bridges.Variable.AbstractBridge
    affine::Vector{MOI.ScalarAffineFunction{T}}
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
    affine = MB.coefficients(acc, set.basis)
    return KernelBridge{T}(affine, variables, constraints)
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
    return [MOI.get(model, MOI.VariablePrimal(attr.result_index), bridge, MOI.Bridges.IndexInVector(i)) for i in eachindex(bridge.affine)]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::KernelBridge,
    i::MOI.Bridges.IndexInVector,
)
    return MOI.Utilities.eval_variable(bridge.affine[i.value]) do
        vi -> MOI.get(model, MOI.VariablePrimal(attr.result_index), vi)
    end
end

function MOI.Bridges.bridged_function(
    bridge::KernelBridge,
    i::MOI.Bridges.IndexInVector,
)
    return bridge.affine[i.value]
end
