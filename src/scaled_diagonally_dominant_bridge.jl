struct ScaledDiagonallyDominantBridge{T, F, VBS} <: MOIB.AbstractBridge
    side_dimension::Int
    variable_bridge::VBS
    equality::MOI.ConstraintIndex{F, MOI.Zeros}
end

function ScaledDiagonallyDominantBridge{T, F, VBS}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::ScaledDiagonallyDominantConeTriangle) where {T, F, VBS}
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    Q, variable_bridge = add_matrix_variable_bridge(
        model, ScaledDiagonallyDominantConeTriangle, side_dimension(s), T)
    g = MOI.operate(-, T, f, MOIU.vectorize(g))
    equality = MOI.add_constraint(model, g, MOI.Zeros(MOI.dimension(s)))
    return ScaledDiagonallyDominantBridge{T, F, VBS}(
        side_dimension(s), variable_bridge, equality)
end

function MOI.supports_constraint(::Type{<:ScaledDiagonallyDominantBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ScaledDiagonallyDominantConeTriangle})
    return true
end
function MOIB.added_constraint_types(::Type{<:ScaledDiagonallyDominantBridge{T}}) where {T}
    added = [(F, MOI.Zeros)]
    return append_added_constraint_types(
        added, ScaledDiagonallyDominantConeTriangle, T)
end
function MOIB.concrete_bridge_type(::Type{<:ScaledDiagonallyDominantBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{ScaledDiagonallyDominantConeTriangle}) where T
    G = MOIU.promote_operation(-, T, F, MOI.VectorAffineFunction{T})
    VBS = union_vector_bridge_types(ScaledDiagonallyDominantConeTriangle, T)
    return ScaledDiagonallyDominantBridge{T, G, VBS}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 attr::MOI.NumberOfVariables)
    return MOI.get(bridge.variable_bridge, attr)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 attr::MOI.NumberOfConstraints)
    return MOI.get(bridge.variable_bridge, attr)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 attr::MOI.ListOfConstraintIndices)
    return MOI.get(bridge.variable_bridge, attr)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.Zeros}) where {T, F}
    return 1
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.Zeros}) where {T, F}
    return [bridge.equality]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::ScaledDiagonallyDominantBridge)
    MOI.delete(model, bridge.equality)
    MOI.delete(model, bridge.variable_bridge)
end

# TODO ConstraintPrimal
function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::ScaledDiagonallyDominantBridge)
    dual = copy(MOI.get(model, attr, bridge.equality))
    # Need to divide by 2 because of the custom scalar product for this cone
    k = 0
    for j in 1:bridge.side_dimension
        for i in 1:(j-1)
            k += 1
            dual[k] /= 2
        end
        k += 1
    end
    return dual
end
