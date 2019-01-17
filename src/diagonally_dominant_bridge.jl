struct DiagonallyDominantBridge{T, F} <: MOIB.AbstractBridge
    # inequalities Qjj ≥ sum_{i ≠ j} |Qij|
    dominance::Vector{MOI.ConstraintIndex{F, MOI.GreaterThan{T}}}
end

function DiagonallyDominantBridge{T, F}(model::MOI.ModelLike,
                                        f::MOI.AbstractVectorFunction,
                                        s::DiagonallyDominantConeTriangle) where T
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    n = s.side_dimension
    g = F[zero(F) for i in 1:n]
    fs = MOIU.eachscalar(f)
    k = 0
    for j in 1:n
        # Qjj ≥ sum_{i ≠ j} |Qij|
        for i in 1:(j-1)
            k += 1
            # abs ≥ |Qij|
            abs = MOI.add_variable(model)
            MOIU.operate!(-, T, dominance[j], abs)
            MOIU.operate!(-, T, dominance[i], abs)
            MOI.add_constraint(model, MOIU.operate(+, T, abs, fs[k]),
                               MOI.GreaterThan(0.0))
            MOI.add_constraint(model, MOIU.operate(-, T, abs, fs[k]),
                               MOI.GreaterThan(0.0))
        end
        k += 1
        MOIU.operate!(+, T, g[j], fs[k])
    end
    dominance = map(f -> MOI.add_constraint(model, f, MOI.GreaterThan(0.0)), g)
    return DiagonallyDominantBridge{T, F}(dominance)
end

function MOI.supports_constraint(::Type{<:DiagonallyDominantBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:DiagonallyDominantConeTriangle})
    return true
end
function MOIB.added_constraint_types(::Type{DiagonallyDominantBridge{T, F}}) where {T, F}
    return [(F, MOI.GreaterThan{T})]
end
function MOIB.concrete_bridge_type(::Type{<:DiagonallyDominantBridge},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{DiagonallyDominantConeTriangle})
    S = MOIU.scalar_type(F)
    G = MOIU.promote_operation(-, T, S, MOI.SingleVariable)
    return DiagonallyDominantBridge{T, G}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::DiagonallyDominantBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.GreaterThan{T}}) where {T, F}
    return length(bridge.dominance)
end
function MOI.get(bridge::DiagonallyDominantBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.GreaterThan{T}}) where {T, F}
    return bridge.dominance
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::DiagonallyDominantBridge)
    for ci in bridge.dominance
        MOI.delete(model, ci)
    end
    # TODO delete variables
end

# TODO ConstraintPrimal and ConstraintDual
