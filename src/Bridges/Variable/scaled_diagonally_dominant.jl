"""
    struct ScaledDiagonallyDominantBridge{T} <: MOI.Bridges.Variable.AbstractBridge
        side_dimension::Int
        variables::Vector{Vector{MOI.VariableIndex}}
        constraints::Vector{MOI.ConstraintIndex{
            MOI.VectorOfVariables, SOS.PositiveSemidefinite2x2ConeTriangle}}
    end

A matrix is SDD iff it is the sum of psd matrices Mij that are zero except
for entries ii, ij and jj [Lemma 9, AM17]. This bridge substitute the
constrained variables in [`SumOfSquares.ScaledDiagonallyDominantConeTriangle`](@ref)
into a sum of constrained variables in [`SumOfSquares.PositiveSemidefinite2x2ConeTriangle`](@ref).

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, **2017**.
"""
struct ScaledDiagonallyDominantBridge{T} <: MOIB.Variable.AbstractBridge
    side_dimension::Int
    variables::Vector{Vector{MOI.VariableIndex}}
    constraints::Vector{MOI.ConstraintIndex{
        MOI.VectorOfVariables, SOS.PositiveSemidefinite2x2ConeTriangle}}
end

function MOIB.Variable.bridge_constrained_variable(
    ::Type{ScaledDiagonallyDominantBridge{T}}, model::MOI.ModelLike,
    s::SOS.ScaledDiagonallyDominantConeTriangle) where T
    n = s.side_dimension
    if n <= 1
        error("The bridges does not work with 1, `matrix_cone` should have returned `Nonnegatives` instead.")
    end
    N = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(n))
    variables = Vector{Vector{MOI.VariableIndex}}(undef, N)
    constraints = Vector{MOI.ConstraintIndex{MOI.VectorOfVariables, SOS.PositiveSemidefinite2x2ConeTriangle}}(undef, N)
    k = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            # PSD constraints on 2x2 matrices are SOC representable
            variables[k], constraints[k] = MOI.add_constrained_variables(model, SOS.PositiveSemidefinite2x2ConeTriangle())
        end
    end
    return ScaledDiagonallyDominantBridge{T}(n, variables, constraints)
end

function MOIB.Variable.supports_constrained_variable(
    ::Type{<:ScaledDiagonallyDominantBridge},
    ::Type{SOS.ScaledDiagonallyDominantConeTriangle})
    return true
end
function MOIB.added_constrained_variable_types(
    ::Type{<:ScaledDiagonallyDominantBridge})
    return [(SOS.PositiveSemidefinite2x2ConeTriangle,)]
end
function MOIB.added_constraint_types(
    ::Type{<:ScaledDiagonallyDominantBridge})
    return Tuple{DataType, DataType}[]
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 ::MOI.NumberOfVariables)
    return 3length(bridge.variables)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 ::MOI.ListOfVariableIndices)
    return Iterators.flatten(bridge.variables)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, SOS.PositiveSemidefinite2x2ConeTriangle})
    return length(bridge.constraints)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge,
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, SOS.PositiveSemidefinite2x2ConeTriangle})
    return bridge.constraints
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::ScaledDiagonallyDominantBridge)
    for variables in bridge.variables
        MOI.delete(model, variables)
    end
end

# Attributes, Bridge acting as a constraint

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet,
                 bridge::ScaledDiagonallyDominantBridge)
    return SOS.ScaledDiagonallyDominantConeTriangle(bridge.side_dimension)
end

# TODO ConstraintPrimal, ConstraintDual

trimap(i, j) = div(j * (j - 1), 2) + i

function MOI.get(model::MOI.ModelLike, attr::MOI.VariablePrimal,
                 bridge::ScaledDiagonallyDominantBridge{T}, i::MOIB.Variable.IndexInVector) where T
    i, j = matrix_indices(i.value)
    if i == j
        value = zero(T)
        for k in 1:(i - 1)
            idx = offdiag_vector_index(k, i)
            value += MOI.get(model, attr, bridge.variables[idx][3])
        end
        for k in (i + 1):bridge.side_dimension
            idx = offdiag_vector_index(i, k)
            value += MOI.get(model, attr, bridge.variables[idx][1])
        end
        return value
    else
        idx = offdiag_vector_index(i, j)
        return MOI.get(model, attr, bridge.variables[idx][2])
    end
end

function MOIB.bridged_function(bridge::ScaledDiagonallyDominantBridge{T},
                               i::MOIB.Variable.IndexInVector) where T
    i, j = matrix_indices(i.value)
    if i == j
        func = zero(MOI.ScalarAffineFunction{T})
        for k in 1:(i - 1)
            idx = offdiag_vector_index(k, i)
            MOIU.operate!(+, T, func, MOI.SingleVariable(bridge.variables[idx][3]))
        end
        for k in (i + 1):bridge.side_dimension
            idx = offdiag_vector_index(i, k)
            MOIU.operate!(+, T, func, MOI.SingleVariable(bridge.variables[idx][1]))
        end
        return func
    else
        idx = offdiag_vector_index(i, j)
        return MOI.convert(MOI.ScalarAffineFunction{T}, MOI.SingleVariable(bridge.variables[idx][2]))
    end
end
function MOIB.Variable.unbridged_map(
    bridge::ScaledDiagonallyDominantBridge{T},
    vis::Vector{MOI.VariableIndex}) where T

    SAF = MOI.ScalarAffineFunction{T}
    umap = Pair{MOI.VariableIndex, SAF}[]
    k = 0
    z = zero(SAF)
    saf(i) = convert(SAF, MOI.SingleVariable(vis[i]))
    # vis[trimap(j, j)] is replaced by a sum of several variables.
    # The strategy is to replace all of them by zero except one.
    for j in 1:bridge.side_dimension
        for i in 1:(j-1)
            k += 1
            if i == 1 && j == 2
                push!(umap, bridge.variables[k][1] => saf(trimap(1, 1)))
            else
                push!(umap, bridge.variables[k][1] => z)
            end
            push!(umap, bridge.variables[k][2] => saf(trimap(i, j)))
            if i == 1
                push!(umap, bridge.variables[k][3] => saf(trimap(j, j)))
            else
                push!(umap, bridge.variables[k][3] => z)
            end
        end
    end
    return umap
end
