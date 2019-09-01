struct DiagonallyDominantBridge{T, F} <: MOIB.Constraint.AbstractBridge
    # |Qij| variables
    abs_vars::Vector{MOI.VariableIndex}
    # |Qij| ≥ +Qij
    abs_plus::Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                         MOI.GreaterThan{T}}}
    # |Qij| ≥ -Qij
    abs_minus::Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},
                                          MOI.GreaterThan{T}}}
    # inequalities Qjj ≥ sum_{i ≠ j} |Qij|
    dominance::Vector{MOI.ConstraintIndex{F, MOI.GreaterThan{T}}}
end

function MOIB.Constraint.bridge_constraint(
    ::Type{DiagonallyDominantBridge{T, F}},
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::DiagonallyDominantConeTriangle) where {T, F}

    @assert MOI.output_dimension(f) == MOI.dimension(s)
    n = s.side_dimension
    g = F[zero(F) for i in 1:n]
    fs = MOIU.eachscalar(f)
    num_off_diag = MOI.dimension(s) - n
    CI = MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, MOI.GreaterThan{T}}
    abs_vars = Vector{MOI.VariableIndex}(undef, num_off_diag)
    abs_plus = Vector{CI}(undef, num_off_diag)
    abs_minus = Vector{CI}(undef, num_off_diag)
    k = 0
    koff = 0
    for j in 1:n
        # Qjj ≥ sum_{i ≠ j} |Qij|
        for i in 1:(j-1)
            k += 1
            koff += 1
            # abs ≥ |Qij|
            abs_vars[koff] = MOI.add_variable(model)
            fabs = MOI.SingleVariable(abs_vars[koff])
            MOIU.operate!(-, T, g[j], fabs)
            MOIU.operate!(-, T, g[i], fabs)
            abs_plus[koff] = MOI.add_constraint(
                model, MOIU.operate(+, T, fabs, fs[k]), MOI.GreaterThan(0.0))
            abs_minus[koff] = MOI.add_constraint(
                model, MOIU.operate(-, T, fabs, fs[k]), MOI.GreaterThan(0.0))
        end
        k += 1
        MOIU.operate!(+, T, g[j], fs[k])
    end
    dominance = map(f -> MOI.add_constraint(model, f, MOI.GreaterThan(0.0)), g)
    return DiagonallyDominantBridge{T, F}(abs_vars, abs_plus, abs_minus,
                                          dominance)
end

function MOI.supports_constraint(::Type{<:DiagonallyDominantBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:DiagonallyDominantConeTriangle})
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:DiagonallyDominantBridge})
    return Tuple{DataType}[]
end
function MOIB.added_constraint_types(::Type{DiagonallyDominantBridge{T, F}}) where {T, F}
    added = [(F, MOI.GreaterThan{T})]
    if F != MOI.ScalarAffineFunction{T}
        push!(added, (MOI.ScalarAffineFunction{T}, MOI.GreaterThan{T}))
    end
    return added
end
function MOIB.Constraint.concrete_bridge_type(
    ::Type{<:DiagonallyDominantBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{DiagonallyDominantConeTriangle}) where T

    S = MOIU.scalar_type(F)
    G = MOIU.promote_operation(-, T, S, MOI.SingleVariable)
    return DiagonallyDominantBridge{T, G}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::DiagonallyDominantBridge, ::MOI.NumberOfVariables)
    return length(bridge.abs_vars)
end
function MOI.get(bridge::DiagonallyDominantBridge{T, MOI.ScalarAffineFunction{T}},
                 ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T},
                                           MOI.GreaterThan{T}}) where T
    return length(bridge.abs_plus) + length(bridge.abs_minus) + length(bridge.dominance)
end
function MOI.get(bridge::DiagonallyDominantBridge{T},
                 ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T},
                                           MOI.GreaterThan{T}}) where T
    return length(bridge.abs_plus) + length(bridge.abs_minus)
end
function MOI.get(bridge::DiagonallyDominantBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.GreaterThan{T}}) where {T, F}
    return length(bridge.dominance)
end
function MOI.get(bridge::DiagonallyDominantBridge{T, MOI.ScalarAffineFunction{T}},
                 ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},
                                               MOI.GreaterThan{T}}) where T
    return vcat(bridge.abs_plus, bridge.abs_minus, bridge.dominance)
end
function MOI.get(bridge::DiagonallyDominantBridge{T},
                 ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},
                                               MOI.GreaterThan{T}}) where T
    return vcat(bridge.abs_plus, bridge.abs_minus)
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
    for ci in bridge.abs_plus
        MOI.delete(model, ci)
    end
    for ci in bridge.abs_minus
        MOI.delete(model, ci)
    end
    for vi in bridge.abs_vars
        MOI.delete(model, vi)
    end
end

# TODO ConstraintPrimal

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::DiagonallyDominantBridge{T}) where T
    dominance_dual = MOI.get(model, attr, bridge.dominance)
    side_dim = length(dominance_dual)
    dim = MOI.dimension(DiagonallyDominantConeTriangle(side_dim))
    dual = Array{T}(undef, dim)
    k = 0
    for j in 1:side_dim
        for i in 1:(j-1)
            k += 1
            # Need to divide by 2 because of the custom scalar product for this
            # cone
            dual[k] = (- dominance_dual[i] - dominance_dual[j]) / 2
        end
        k += 1
        dual[k] = dominance_dual[j]
    end
    @assert k == dim
    return dual
end
