"""
    struct ScaledDiagonallyDominantBridge{T} <: MOI.Bridges.Variable.AbstractBridge
        side_dimension::Int
        variables::Vector{Vector{MOI.VariableIndex}}
        constraints::Vector{MOI.ConstraintIndex{
            MOI.VectorOfVariables, SOS.PositiveSemidefinite2x2ConeTriangle}}
    end

A matrix is SDD iff it is the sum of psd matrices Mij that are zero except
for entries ii, ij and jj [Ahmadi2017; Lemma 9](@cite). This bridge substitute the
constrained variables in [`SOS.ScaledDiagonallyDominantConeTriangle`](@ref)
into a sum of constrained variables in [`SOS.PositiveSemidefinite2x2ConeTriangle`](@ref).
"""
struct ScaledDiagonallyDominantBridge{T} <: MOI.Bridges.Variable.AbstractBridge
    side_dimension::Int
    variables::Vector{Vector{MOI.VariableIndex}}
    constraints::Vector{
        MOI.ConstraintIndex{
            MOI.VectorOfVariables,
            SOS.PositiveSemidefinite2x2ConeTriangle,
        },
    }
end

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{ScaledDiagonallyDominantBridge{T}},
    model::MOI.ModelLike,
    s::SOS.ScaledDiagonallyDominantConeTriangle,
) where {T}
    n = s.side_dimension
    if n <= 1
        error(
            "The bridges does not work with 1, `matrix_cone` should have returned `Nonnegatives` instead.",
        )
    end
    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n - 1))
    variables = Vector{Vector{MOI.VariableIndex}}(undef, N)
    constraints = Vector{
        MOI.ConstraintIndex{
            MOI.VectorOfVariables,
            SOS.PositiveSemidefinite2x2ConeTriangle,
        },
    }(
        undef,
        N,
    )
    k = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            # PSD constraints on 2x2 matrices are SOC representable
            variables[k], constraints[k] = MOI.add_constrained_variables(
                model,
                SOS.PositiveSemidefinite2x2ConeTriangle(),
            )
        end
    end
    return ScaledDiagonallyDominantBridge{T}(n, variables, constraints)
end

function MOI.Bridges.Variable.supports_constrained_variable(
    ::Type{<:ScaledDiagonallyDominantBridge},
    ::Type{SOS.ScaledDiagonallyDominantConeTriangle},
)
    return true
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ScaledDiagonallyDominantBridge},
)
    return Tuple{Type}[(SOS.PositiveSemidefinite2x2ConeTriangle,)]
end
function MOI.Bridges.added_constraint_types(
    ::Type{<:ScaledDiagonallyDominantBridge},
)
    return Tuple{Type,Type}[]
end

# Attributes, Bridge acting as a model
function MOI.get(
    bridge::ScaledDiagonallyDominantBridge,
    ::MOI.NumberOfVariables,
)
    return 3length(bridge.variables)
end
function MOI.get(
    bridge::ScaledDiagonallyDominantBridge,
    ::MOI.ListOfVariableIndices,
)
    return collect(Iterators.flatten(bridge.variables))
end
function MOI.get(
    bridge::ScaledDiagonallyDominantBridge,
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        SOS.PositiveSemidefinite2x2ConeTriangle,
    },
)
    return length(bridge.constraints)
end
function MOI.get(
    bridge::ScaledDiagonallyDominantBridge,
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        SOS.PositiveSemidefinite2x2ConeTriangle,
    },
)
    return bridge.constraints
end

# Indices
function MOI.delete(
    model::MOI.ModelLike,
    bridge::ScaledDiagonallyDominantBridge,
)
    for variables in bridge.variables
        MOI.delete(model, variables)
    end
end

# Attributes, Bridge acting as a constraint

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::ScaledDiagonallyDominantBridge,
)
    return SOS.ScaledDiagonallyDominantConeTriangle(bridge.side_dimension)
end

# The map `A` is not injective because it maps the entry of several 2x2 matrices
# into the same index so it is not invertible hence it's unclear how to implemented
# `set` for `VariablePrimalStart`.
# The adjoint `A'` is however injective and it is easy to invert.

function MOI.supports(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDualStart,
    ::Type{<:ScaledDiagonallyDominantBridge},
)
    return MOI.supports(
        model,
        attr,
        MOI.ConstraintIndex{
            MOI.VectorOfVariables,
            SOS.PositiveSemidefinite2x2ConeTriangle,
        },
    )
end

function MOI.set(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDualStart,
    bridge::ScaledDiagonallyDominantBridge,
    value,
)
    n = bridge.side_dimension
    k = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            # PSD constraints on 2x2 matrices are SOC representable
            if isnothing(value)
                dual = nothing
            else
                dual = [
                    value[MOI.Utilities.trimap(i, i)],
                    value[MOI.Utilities.trimap(i, j)],
                    value[MOI.Utilities.trimap(j, j)],
                ]
            end
            MOI.set(model, attr, bridge.constraints[k], dual)
        end
    end
    return
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{MOI.ConstraintDual,MOI.ConstraintDualStart},
    bridge::ScaledDiagonallyDominantBridge{T},
) where {T}
    n = bridge.side_dimension
    value = zeros(T, MOI.Utilities.trimap(n, n))
    k = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            dual = MOI.get(model, attr, bridge.constraints[k])
            if isnothing(dual)
                return nothing
            end
            # There are `bridge.side_dimension - 1` possible candidate that should all have
            # the same `dual` so we take an arbitrary choice
            if j == i + 1
                value[MOI.Utilities.trimap(i, i)] = dual[1]
            elseif i == 1 && j == n
                value[MOI.Utilities.trimap(j, j)] = dual[3]
            end
            value[MOI.Utilities.trimap(i, j)] = dual[2]
        end
    end
    return value
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::ScaledDiagonallyDominantBridge{T},
    index::MOI.Bridges.IndexInVector,
) where {T}
    i, j = MOI.Utilities.inverse_trimap(index.value)
    if i == j
        value = zero(T)
        for k in 1:(i-1)
            idx = offdiag_vector_index(k, i)
            value += MOI.get(model, attr, bridge.variables[idx][3])
        end
        for k in (i+1):bridge.side_dimension
            idx = offdiag_vector_index(i, k)
            value += MOI.get(model, attr, bridge.variables[idx][1])
        end
        return value
    else
        idx = offdiag_vector_index(i, j)
        return MOI.get(model, attr, bridge.variables[idx][2])
    end
end

function MOI.Bridges.bridged_function(
    bridge::ScaledDiagonallyDominantBridge{T},
    i::MOI.Bridges.IndexInVector,
) where {T}
    i, j = MOI.Utilities.inverse_trimap(i.value)
    if i == j
        func = zero(MOI.ScalarAffineFunction{T})
        for k in 1:(i-1)
            idx = offdiag_vector_index(k, i)
            MOI.Utilities.operate!(+, T, func, bridge.variables[idx][3])
        end
        for k in (i+1):bridge.side_dimension
            idx = offdiag_vector_index(i, k)
            MOI.Utilities.operate!(+, T, func, bridge.variables[idx][1])
        end
        return func
    else
        idx = offdiag_vector_index(i, j)
        return MOI.convert(
            MOI.ScalarAffineFunction{T},
            bridge.variables[idx][2],
        )
    end
end
function MOI.Bridges.Variable.unbridged_map(
    bridge::ScaledDiagonallyDominantBridge{T},
    vis::Vector{MOI.VariableIndex},
) where {T}
    SAF = MOI.ScalarAffineFunction{T}
    umap = Pair{MOI.VariableIndex,SAF}[]
    k = 0
    z = zero(SAF)
    saf(i) = convert(SAF, vis[i])
    # vis[MOI.Utilities.trimap(j, j)] is replaced by a sum of several variables.
    # The strategy is to replace all of them by zero except one.
    for j in 1:bridge.side_dimension
        for i in 1:(j-1)
            k += 1
            if i == 1 && j == 2
                push!(
                    umap,
                    bridge.variables[k][1] => saf(MOI.Utilities.trimap(1, 1)),
                )
            else
                push!(umap, bridge.variables[k][1] => z)
            end
            push!(
                umap,
                bridge.variables[k][2] => saf(MOI.Utilities.trimap(i, j)),
            )
            if i == 1
                push!(
                    umap,
                    bridge.variables[k][3] => saf(MOI.Utilities.trimap(j, j)),
                )
            else
                push!(umap, bridge.variables[k][3] => z)
            end
        end
    end
    return umap
end
