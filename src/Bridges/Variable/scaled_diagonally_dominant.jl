struct ScaledDiagonallyDominantVariableBridge{T, VBT} <: AbstractVariableBridge
    psd2x2::Vector{VBT}
end

function add_variable_bridge(
    ::Type{ScaledDiagonallyDominantVariableBridge{T, VBT}}, model::MOI.ModelLike,
    s::ScaledDiagonallyDominantConeTriangle) where {T, VBT}
    # `p.Q` is SDD iff it is the sum of psd matrices Mij that are zero except
    # for entries ii, ij and jj [Lemma 9, AM17].
    n = s.side_dimension
    @assert n > 1 # The bridges does not work with 1, `matrix_cone`
                  # should have returned `Nonnegatives` instead
    # `g[r, c]` will contain the expression `sum Mij[r, c]`.
    F = MOI.ScalarAffineFunction{T}
    g = F[zero(F) for i in 1:MOI.dimension(s)]
    num_off_diag = MOI.dimension(s) - n
    psd2x2 = Vector{VBT}(undef, num_off_diag)
    diag_idx(i) = div(i * (i + 1), 2)
    k = 0
    k2x2 = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            k2x2 += 1
            # PSD constraints on 2x2 matrices are SOC representable
            x, psd2x2[k2x2] = add_variable_bridge(VBT, model, PositiveSemidefinite2x2ConeTriangle())
            MOIU.operate!(+, T, g[diag_idx(i)], x[1])
            MOIU.operate!(+, T, g[k], x[2])
            MOIU.operate!(+, T, g[diag_idx(j)], x[3])
        end
        k += 1 # diagonal entry `(i, i)`
    end
    return g, ScaledDiagonallyDominantVariableBridge{T, VBT}(psd2x2)
end

function MOIB.added_constraint_types(
    ::Type{ScaledDiagonallyDominantVariableBridge{T, VBT}}) where {T, VBT}
    added = Tuple{DataType, DataType}[]
    return append!(added, MOIB.added_constraint_types(VBT))
end

function variable_bridge_type(::Type{ScaledDiagonallyDominantConeTriangle},
                              T::Type)
    VBT = variable_bridge_type(PositiveSemidefinite2x2ConeTriangle, T)
    return ScaledDiagonallyDominantVariableBridge{T, VBT}
end


# Attributes, VariableBridge acting as an model
function MOI.get(bridge::ScaledDiagonallyDominantVariableBridge, attr::MOI.NumberOfVariables)
    return reduce(+, MOI.get.(bridge.psd2x2, attr), init=0)
end
function MOI.get(bridge::ScaledDiagonallyDominantVariableBridge,
                 attr::MOI.NumberOfConstraints)
    return reduce(+, MOI.get.(bridge.psd2x2, attr), init=0)
end
function MOI.get(bridge::ScaledDiagonallyDominantVariableBridge,
                 attr::MOI.ListOfConstraintIndices{F, S}) where {F, S}
    list = MOI.ConstraintIndex{F, S}[]
    for variable_bridge in bridge.psd2x2
        append!(list, MOI.get(variable_bridge, attr))
    end
    return list
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::ScaledDiagonallyDominantVariableBridge)
    for ci in bridge.psd2x2
        MOI.delete(model, ci)
    end
end
