struct LowRankBridge{T,M,B,G,W} <: MOI.Bridges.Variable.AbstractBridge
    affine::Vector{MOI.ScalarAffineFunction{T}}
    variables::Vector{Vector{MOI.VariableIndex}}
    constraints::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables}}
    set::SOS.WeightedSOSCone{M,B,G,W}
end

import LinearAlgebra

# The matrix returned by `MB.transformation_to(gram_basis, target_basis)`.
# `Base.promote_op` resolves this at the type level so we can declare a
# concrete `LRO.Factorization{T, SubArray{...}, Array{T,0}}` in
# `added_constrained_variable_types` — needed because PolyJuMP's
# `bridgeable` rejects `UnionAll` set types.
_transformation_type(::Type{G}, ::Type{B}) where {G,B} =
    Base.promote_op(MB.transformation_to, G, B)

# Row-view of the transformation matrix, i.e. the type of `view(U, j, :)`.
# Used as the `Factorization.factor` type so the parent `U` survives down
# to downstream batched-FFT consumers (e.g. `LowRankOpt.BurerMonteiro`).
_row_view_type(::Type{MT}) where {MT} =
    Base.promote_op(view, MT, Int, Colon)

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{LowRankBridge{T,M,B,G,W}},
    model::MOI.ModelLike,
    set::SOS.WeightedSOSCone{M,B,G,W},
) where {T,M,B,G,W}
    variables = Vector{Vector{MOI.VariableIndex}}(undef, length(set.gram_bases))
    constraints = Vector{MOI.ConstraintIndex{MOI.VectorOfVariables}}(
        undef,
        length(set.gram_bases),
    )
    for i in eachindex(set.gram_bases)
        U = MB.transformation_to(set.gram_bases[i], set.basis)
        weights = SA.coeffs(set.weights[i], set.basis)
        # `view(U, j, :)` preserves the parent matrix `U` so that downstream
        # `LowRankOpt.BurerMonteiro` can recognise the row-sharing pattern and
        # use a single `mul!(buf, U, X.factor)` (one batched FFT per row of `Y`
        # when `U::MB.TrigEvalMatrix`) instead of `n` separate inner products.
        variables[i], constraints[i] = MOI.add_constrained_variables(
            model,
            LRO.SetDotProducts{LRO.WITHOUT_SET}(
                SOS.matrix_cone(M, length(set.gram_bases[i])),
                [
                    LRO.TriangleVectorization(
                        LRO.Factorization(
                            view(U, j, :),
                            reshape(T[weights[j]], ()),
                        ),
                    ) for j in eachindex(set.basis)
                ],
            ),
        )
    end
    return LowRankBridge{T,M,B,G,W}(
        [
            MOI.ScalarAffineFunction(
                [
                    MOI.ScalarAffineTerm(one(T), variables[i][j]) for
                    i in eachindex(set.gram_bases)
                ],
                zero(T),
            ) for j in eachindex(set.basis)
        ],
        variables,
        constraints,
        set,
    )
end

function MOI.Bridges.Variable.supports_constrained_variable(
    ::Type{<:LowRankBridge},
    ::Type{<:SOS.WeightedSOSCone{M,B}},
) where {M,B}
    # Could be made to work for non-LagrangeBasis but it's not low rank in
    # we we'll need a high bridge cost for them so that it's only UnsafeAddMul
    # if the other bridges are removed
    return B <: MB.LagrangeBasis
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{LowRankBridge{T,M,B,G,W}},
) where {T,M,B,G,W}
    MT = _transformation_type(G, B)
    FT = _row_view_type(MT)
    TVT = LRO.TriangleVectorization{T,LRO.Factorization{T,FT,Array{T,0}}}
    return Tuple{Type}[
        (LRO.SetDotProducts{LRO.WITHOUT_SET,S[1],TVT,Vector{TVT}},)
        for S in SOS.Bridges.Constraint.constrained_variable_types(M) if
        S[1] == MOI.PositiveSemidefiniteConeTriangle # FIXME hack
    ]
end

function MOI.Bridges.added_constraint_types(::Type{<:LowRankBridge})
    return Tuple{Type,Type}[]
end

function MOI.Bridges.Variable.concrete_bridge_type(
    ::Type{<:LowRankBridge{T}},
    ::Type{<:SOS.WeightedSOSCone{M,B,G,W}},
) where {T,M,B,G,W}
    return LowRankBridge{T,M,B,G,W}
end

# Attributes, Bridge acting as a model
function MOI.get(bridge::LowRankBridge, ::MOI.NumberOfVariables)
    return sum(length, bridge.variables)
end

function MOI.get(bridge::LowRankBridge, ::MOI.ListOfVariableIndices)
    return reduce(vcat, bridge.variables)
end

function MOI.get(
    bridge::LowRankBridge,
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
) where {S<:MOI.AbstractVectorSet}
    return count(bridge.constraints) do ci
        return ci isa MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    end
end

function MOI.get(
    bridge::LowRankBridge,
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {S}
    return [
        ci for ci in bridge.constraints if
        ci isa MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    ]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::LowRankBridge)
    for vars in bridge.variables
        MOI.delete(model, vars)
    end
    return
end

# Attributes, Bridge acting as a constraint

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, bridge::LowRankBridge)
    return bridge.set
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintPrimal,
    bridge::LowRankBridge,
)
    return [
        MOI.get(
            model,
            MOI.VariablePrimal(attr.result_index),
            bridge,
            MOI.Bridges.IndexInVector(i),
        ) for i in eachindex(bridge.affine)
    ]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.VariablePrimal,
    bridge::LowRankBridge,
    i::MOI.Bridges.IndexInVector,
)
    return MOI.Utilities.eval_variables(bridge.affine[i.value]) do vi
        return MOI.get(model, MOI.VariablePrimal(attr.result_index), vi)
    end
end

function MOI.Bridges.bridged_function(
    bridge::LowRankBridge,
    i::MOI.Bridges.IndexInVector,
)
    return bridge.affine[i.value]
end

function MOI.Bridges.Variable.unbridged_map(
    ::LowRankBridge,
    ::Vector{MOI.VariableIndex},
)
    return nothing
end
