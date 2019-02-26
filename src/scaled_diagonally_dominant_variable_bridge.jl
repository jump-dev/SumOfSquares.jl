struct ScaledDiagonallyDominantVariableBridge{T} <: AbstractVariableBridge
    side_dimension::Int
    variables::Vector{NTuple{3, MOI.VariableIndex}}
    psd2x2::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables,
                                       PositiveSemidefinite2x2ConeTriangle}}
end

function add_variable_bridge(
    ::Type{ScaledDiagonallyDominantVariableBridge{T}}, model::MOI.ModelLike,
    s::ScaledDiagonallyDominantConeTriangle) where {T}
    # `p.Q` is SDD iff it is the sum of psd matrices Mij that are zero except
    # for entries ii, ij and jj [Lemma 9, AM17].
    n = s.side_dimension
    @assert n > 1 # The bridges does not work with 1, `matrix_cone`
                  # should have returned `Nonnegatives` instead
    # `g[r, c]` will contain the expression `sum Mij[r, c]`.
    F = MOI.ScalarAffineFunction{T}
    g = F[zero(F) for i in 1:MOI.dimension(s)]
    num_off_diag = MOI.dimension(s) - n
    variables = Vector{NTuple{3, MOI.VariableIndex}}(undef, num_off_diag)
    psd2x2 = Vector{MOI.ConstraintIndex{
        MOI.VectorOfVariables, PositiveSemidefinite2x2ConeTriangle}}(
            undef, num_off_diag)
    diag_idx(i) = div(i * (i + 1), 2)
    k = 0
    k2x2 = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            k2x2 += 1
            vii, vij, vjj = MOI.add_variables(model, 3)
            variables[k2x2] = (vii, vij, vjj)
            Mii = MOI.SingleVariable(vii)
            MOIU.operate!(+, T, g[diag_idx(i)], Mii)
            Mij = MOI.SingleVariable(vij)
            MOIU.operate!(+, T, g[k], Mij)
            Mjj = MOI.SingleVariable(vjj)
            MOIU.operate!(+, T, g[diag_idx(j)], Mjj)
            # PSD constraints on 2x2 matrices are SOC representable
            psd2x2[k2x2] = MOI.add_constraint(
                model, MOI.VectorOfVariables([vii, vij, vjj]),
                PositiveSemidefinite2x2ConeTriangle())
        end
        k += 1 # diagonal entry `(i, i)`
    end
    return g, ScaledDiagonallyDominantVariableBridge{T}(n, variables, psd2x2)
end

function MOIB.added_constraint_types(
    ::Type{ScaledDiagonallyDominantVariableBridge{T}}) where {T}
    return [(MOI.VectorOfVariables, PositiveSemidefinite2x2ConeTriangle)]
end

function variable_bridge_type(::Type{ScaledDiagonallyDominantConeTriangle},
                              T::Type)
    return ScaledDiagonallyDominantVariableBridge{T}
end


# Attributes, VariableBridge acting as an model
function MOI.get(bridge::ScaledDiagonallyDominantVariableBridge,
                 ::MOI.NumberOfVariables)
    return 3 * length(bridge.variables)
end
function MOI.get(
    bridge::ScaledDiagonallyDominantVariableBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,
                              ScaledDiagonallyDominantConeTriangle})
    return length(bridge.psd2x2)
end
function MOI.get(
    bridge::ScaledDiagonallyDominantVariableBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,
                                  ScaledDiagonallyDominantConeTriangle})
    return bridge.psd2x2
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::ScaledDiagonallyDominantVariableBridge)
    for ci in bridge.psd2x2
        MOI.delete(model, ci)
    end
    for (vii, vij, vjj) in bridge.variables
        MOI.delete(model, vii)
        MOI.delete(model, vij)
        MOI.delete(model, vjj)
    end
end
