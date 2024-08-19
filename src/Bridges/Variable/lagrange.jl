struct LagrangeBridge{T,M} <: MOI.Bridges.Variable.AbstractBridge
    affine::Vector{MOI.ScalarAffineFunction{T}}
    variables::Vector{Vector{MOI.VariableIndex}}
    constraints::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables}}
    set::SOS.WeightedSOSCone{M}
end

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{LagrangeBridge{T,M}},
    model::MOI.ModelLike,
    set::SOS.WeightedSOSCone{M},
) where {T,M}
    variables = Vector{MOI.VariableIndex}[]
    constraints = MOI.ConstraintIndex{MOI.VectorOfVariables}[]
    acc = zero(
        MOI.ScalarAffineFunction{T},
        MB.algebra(MB.implicit_basis(set.basis)),
    )
    for (gram_basis, weight) in zip(set.gram_bases, set.weights)
        gram, vars, con = SOS.add_gram_matrix(model, M, gram_basis, T)
        push!(variables, vars)
        push!(constraints, con)
        MA.operate!(SA.UnsafeAddMul(*), acc, gram, weight)
    end
    MA.operate!(SA.canonical, SA.coeffs(acc))
    return LagrangeBridge{T,M}(
        SA.coeffs(acc, set.basis),
        variables,
        constraints,
        set,
    )
end

function MOI.Bridges.Variable.supports_constrained_variable(
    ::Type{<:LagrangeBridge},
    ::Type{<:SOS.WeightedSOSCone},
)
    return true
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{LagrangeBridge{T,M}},
) where {T,M}
    return SOS.Bridges.Constraint.constrained_variable_types(M)
end

function MOI.Bridges.added_constraint_types(::Type{<:LagrangeBridge})
    return Tuple{Type,Type}[]
end

function MOI.Bridges.Variable.concrete_bridge_type(
    ::Type{<:LagrangeBridge{T}},
    ::Type{<:SOS.WeightedSOSCone{M}},
) where {T,M}
    return LagrangeBridge{T,M}
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::LagrangeBridge, ::MOI.NumberOfVariables)
    return sum(length, bridge.variables)
end

function MOI.get(bridge::LagrangeBridge, ::MOI.ListOfVariableIndices)
    return reduce(vcat, bridge.variables)
end

function MOI.get(
    bridge::LagrangeBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
) where {S<:MOI.AbstractVectorSet}
    return count(bridge.constraints) do ci
        return ci isa MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    end
end

function MOI.get(
    bridge::LagrangeBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {S}
    return [
        ci for ci in bridge.constraints if
        ci isa MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    ]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::LagrangeBridge)
    for vars in bridge.variables
        MOI.delete(model, vars)
    end
    return
end

# Attributes, Bridge acting as a constraint

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, bridge::LagrangeBridge)
    return bridge.set
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintPrimal,
    bridge::LagrangeBridge,
)
    return [
        MOI.get(
            model,
            MOI.VariablePrimal(attr.result_index),
            bridge,
            MOI.Bridges.IndexInVector(i),
        ) for i in eachindex(bridge.affine)
    ]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::LagrangeBridge,
    i::MOI.Bridges.IndexInVector,
)
    return MOI.Utilities.eval_variables(bridge.affine[i.value]) do vi
        return MOI.get(model, MOI.VariablePrimal(attr.result_index), vi)
    end
end

function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.GramMatrixAttribute,
    bridge::LagrangeBridge{T,M},
) where {T,M}
    SOS.check_multiplier_index_bounds(attr, eachindex(bridge.constraints))
    return SOS.build_gram_matrix(
        convert(
            Vector{T},
            MOI.get(
                model,
                MOI.VariablePrimal(attr.result_index),
                bridge.variables[attr.multiplier_index],
            ),
        ),
        bridge.set.gram_bases[attr.multiplier_index],
        M,
        T,
    )
end

function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.MomentMatrixAttribute,
    bridge::LagrangeBridge{T,M},
) where {T,M}
    SOS.check_multiplier_index_bounds(attr, eachindex(bridge.constraints))
    return SOS.build_moment_matrix(
        convert(
            Vector{T},
            MOI.get(
                model,
                MOI.ConstraintDual(attr.result_index),
                bridge.constraints[attr.multiplier_index],
            ),
        ),
        bridge.set.gram_bases[attr.multiplier_index],
    )
end

function MOI.Bridges.bridged_function(
    bridge::LagrangeBridge,
    i::MOI.Bridges.IndexInVector,
)
    return bridge.affine[i.value]
end

function MOI.Bridges.Variable.unbridged_map(
    ::LagrangeBridge{T},
    ::Vector{MOI.VariableIndex},
) where {T}
    return nothing
end
