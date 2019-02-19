struct SOSLikeVariableBridge{T, S} <: AbstractVariableBridge
    gram_matrix::Vector{MOI.VariableIndex}
    gram_constraint::Union{MOI.ConstraintIndex{MOI.VectorOfVariables, S},
                           MOI.ConstraintIndex{MOI.VectorOfVariables, MOI.Nonnegatives},
                           MOI.ConstraintIndex{MOI.VectorOfVariables, PositiveSemidefinite2x2ConeTriangle}}
end

function add_variable_bridge(::Type{SOSLikeVariableBridge{T, S}},
                             model::MOI.ModelLike, set::S) where {T, S}
    gram_matrix = MOI.add_variables(model, MOI.dimension(set))
    func = MOI.SingleVariable[MOI.SingleVariable(v) for v in gram_matrix]
    ci = matrix_add_constraint(model, func, set)
    # Need to specify `S` as it cannot be inferred if, e.g., `ci` is a
    # `MOI.Nonnegatives` constraint
    return func, SOSLikeVariableBridge{T, S}(gram_matrix, ci)
end

function MOIB.added_constraint_types(::Type{SOSLikeVariableBridge{T, S}}) where {T, S}
    added = [(MOI.VectorOfVariables, S),
             (MOI.VectorOfVariables, MOI.Nonnegatives)]
    if S != DiagonallyDominantConeTriangle
        push!(added, (MOI.VectorOfVariables,
                      PositiveSemidefinite2x2ConeTriangle))
    end
    return added
end

function variable_bridge_type(S::Type{<:MOI.AbstractVectorSet}, T::Type)
    return SOSLikeVariableBridge{T, S}
end

function MOI.delete(model::MOI.ModelLike, bridge::SOSLikeVariableBridge)
    # First delete the constraint in which the Gram matrix appears
    MOI.delete(model, bridge.gram_constraint)
    # Now we delete the Gram matrix
    for variable in bridge.gram_matrix
        MOI.delete(model, variable)
    end
end

function MOI.get(model::MOI.ModelLike, ::MomentMatrixAttribute,
                 bridge::SOSLikeVariableBridge)
    return MOI.get(model, MOI.ConstraintDual(), bridge.gram_constraint)
end

function MOI.get(model::MOI.ModelLike, ::GramMatrixAttribute,
                 bridge::SOSLikeVariableBridge)
    return MOI.get(model, MOI.VariablePrimal(), bridge.gram_matrix)
end
