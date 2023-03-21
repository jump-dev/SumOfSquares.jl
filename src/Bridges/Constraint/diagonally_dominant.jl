struct DiagonallyDominantBridge{T, F, G} <: MOI.Bridges.Constraint.AbstractBridge
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

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{DiagonallyDominantBridge{T, F, G}},
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::SOS.DiagonallyDominantConeTriangle) where {T, F, G}

    @assert MOI.output_dimension(f) == MOI.dimension(s)
    n = s.side_dimension
    g = F[zero(F) for i in 1:n]
    fs = MOI.Utilities.eachscalar(f)
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
            fabs = abs_vars[koff]
            MOI.Utilities.operate!(-, T, g[j], fabs)
            MOI.Utilities.operate!(-, T, g[i], fabs)
            abs_plus[koff] = MOI.add_constraint(
                model, MOI.Utilities.operate(+, T, fabs, fs[k]), MOI.GreaterThan(0.0))
            abs_minus[koff] = MOI.add_constraint(
                model, MOI.Utilities.operate(-, T, fabs, fs[k]), MOI.GreaterThan(0.0))
        end
        k += 1
        MOI.Utilities.operate!(+, T, g[j], fs[k])
    end
    dominance = map(f -> MOI.add_constraint(model, f, MOI.GreaterThan(0.0)), g)
    return DiagonallyDominantBridge{T, F, G}(abs_vars, abs_plus, abs_minus,
                                          dominance)
end

function MOI.supports_constraint(::Type{<:DiagonallyDominantBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOS.DiagonallyDominantConeTriangle})
    return true
end
function MOI.Bridges.added_constrained_variable_types(::Type{<:DiagonallyDominantBridge})
    return Tuple{DataType}[]
end
function MOI.Bridges.added_constraint_types(::Type{<:DiagonallyDominantBridge{T, F}}) where {T, F}
    added = [(F, MOI.GreaterThan{T})]
    if F != MOI.ScalarAffineFunction{T}
        push!(added, (MOI.ScalarAffineFunction{T}, MOI.GreaterThan{T}))
    end
    return added
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:DiagonallyDominantBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.DiagonallyDominantConeTriangle}) where T

    S = MOI.Utilities.scalar_type(G)
    F = MOI.Utilities.promote_operation(-, T, S, MOI.VariableIndex)
    return DiagonallyDominantBridge{T, F, G}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::DiagonallyDominantBridge, ::MOI.NumberOfVariables)
    return length(bridge.abs_vars)
end
function MOI.get(bridge::DiagonallyDominantBridge, ::MOI.ListOfVariableIndices)
    return bridge.abs_vars
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

# Attributes, Bridge acting as a constraint
function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet,
                 bridge::DiagonallyDominantBridge)
    return SOS.DiagonallyDominantConeTriangle(length(bridge.dominance))
end
function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintFunction,
                 bridge::DiagonallyDominantBridge{T, F, G}) where {T, F, G}
    set = MOI.get(model, MOI.ConstraintSet(), bridge)
    H = MOI.Utilities.scalar_type(G)
    g = Vector{H}(undef, MOI.dimension(set))
    k = 0
    koff = 0
    for j in 1:MOI.side_dimension(set)
        for i in 1:(j - 1)
            k += 1
            koff += 1
            func = MOI.get(model, attr, bridge.abs_plus[koff])
            g[k] = MOI.Utilities.convert_approx(H, MOI.Utilities.remove_variable(
                func, bridge.abs_vars))
        end
        k += 1
        func = MOI.get(model, attr, bridge.dominance[j])
        g[k] = MOI.Utilities.convert_approx(H, MOI.Utilities.remove_variable(
            func, bridge.abs_vars))
    end
    return MOI.Utilities.vectorize(g)
end

# TODO ConstraintPrimal

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::DiagonallyDominantBridge{T}) where T
    dominance_dual = MOI.get(model, attr, bridge.dominance)
    side_dim = length(dominance_dual)
    dim = MOI.dimension(SOS.DiagonallyDominantConeTriangle(side_dim))
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
