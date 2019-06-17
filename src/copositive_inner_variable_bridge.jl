struct CopositiveInnerVariableBridge{T, VB} <: AbstractVariableBridge
    variable_bridge::VB
    # TODO store GreaterThan constraints
end

function add_variable_bridge(::Type{CopositiveInnerVariableBridge{T, VB}},
                             model::MOI.ModelLike, set::CopositiveInner) where {T, VB}
    func, variable_bridge = add_variable_bridge(VB, model, set.psd_inner)
    F = MOI.ScalarAffineFunction{T}
    # If `func isa F` then `convert` won't create a copy so we explictly copy it.
    # It is less costly to copy it before, e.g. if `q is MOI.SingleVariable`.
    Q = F[convert(F, copy(q)) for q in func]
    k = 0
    for j in 1:side_dimension(set)
        for i in 1:(j-1)
            k += 1
            # TODO: these should be deleted when the bridge is deleted
            Nij = MOI.SingleVariable(MOI.add_variable(model))
            MOI.add_constraint(model, Nij, MOI.GreaterThan(zero(T)))
            MOIU.operate!(+, T, Q[k], Nij)
        end
        k += 1 # for diagonal entry (i, i)
    end
    return Q, CopositiveInnerVariableBridge{T, VB}(variable_bridge)
end

function MOIB.added_constraint_types(::Type{CopositiveInnerVariableBridge{T, VB}}) where {T, VB}
    added = MOIB.added_constraint_types(VB)
    push!(added, (MOI.SingleVariable, MOI.GreaterThan{T}))
    return added
end

function variable_bridge_type(::Type{CopositiveInner{S}}, T::Type) where S
    return CopositiveInnerVariableBridge{T, variable_bridge_type(S, T)}
end

function MOI.delete(model::MOI.ModelLike, bridge::CopositiveInnerVariableBridge)
    # TODO remove GreaterThan constraints
    MOI.delete(model, bridge.variable_bridge)
end

function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintDual, MOI.ConstraintPrimal},
                 bridge::CopositiveInnerVariableBridge)
    return MOI.get(model, attr, bridge.variable_bridge)
end
