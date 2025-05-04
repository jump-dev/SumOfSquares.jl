struct LowRankBridge{T,M} <: MOI.Bridges.Variable.AbstractBridge
    affine::Vector{MOI.ScalarAffineFunction{T}}
    variables::Vector{Vector{MOI.VariableIndex}}
    constraints::Vector{MOI.ConstraintIndex{MOI.VectorOfVariables}}
    set::SOS.WeightedSOSCone{M}
end

import LinearAlgebra

function MOI.Bridges.Variable.bridge_constrained_variable(
    ::Type{LowRankBridge{T,M}},
    model::MOI.ModelLike,
    set::SOS.WeightedSOSCone{M},
) where {T,M}
    variables = Vector{Vector{MOI.VariableIndex}}(undef, length(set.gram_bases))
    constraints = Vector{MOI.ConstraintIndex{MOI.VectorOfVariables}}(undef, length(set.gram_bases))
    for i in eachindex(set.gram_bases)
        U = MB.transformation_to(set.gram_bases[i], set.basis)
        weights = SA.coeffs(set.weights[i], set.basis)
        variables[i], constraints[i] = MOI.add_constrained_variables(
            model,
            MOI.SetWithDotProducts(
                SOS.matrix_cone(M, length(set.gram_bases[i])),
                [
                    MOI.TriangleVectorization(
                        MOI.LowRankMatrix(
                            [weights[j]],
                            reshape(U[j, :], size(U, 2), 1),
                        )
                    )
                    for j in eachindex(set.basis)
                ],
            ),
        )
    end
    return LowRankBridge{T,M}(
        [
            MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(one(T), variables[i][j]) for i in eachindex(set.gram_bases)],
                zero(T),
            )
            for j in eachindex(set.basis)
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
    ::Type{LowRankBridge{T,M}},
) where {T,M}
    return Tuple{Type}[
        (MOI.SetWithDotProducts{S[1],MOI.TriangleVectorization{MOI.LowRankMatrix{T}}},)
        for S in SOS.Bridges.Constraint.constrained_variable_types(M)
        if S[1] == MOI.PositiveSemidefiniteConeTriangle # FIXME hack
    ]
end

function MOI.Bridges.added_constraint_types(::Type{<:LowRankBridge})
    return Tuple{Type,Type}[]
end

function MOI.Bridges.Variable.concrete_bridge_type(
    ::Type{<:LowRankBridge{T}},
    ::Type{<:SOS.WeightedSOSCone{M}},
) where {T,M}
    return LowRankBridge{T,M}
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
