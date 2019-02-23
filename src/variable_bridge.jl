abstract type AbstractVariableBridge end
MOI.get(::AbstractVariableBridge, ::MOI.NumberOfConstraints) = 0
function MOI.get(::AbstractVariableBridge,
                 ::MOI.ListOfConstraintIndices{F, S}) where {F, S}
    return MOI.ConstraintIndex{F, S}[]
end

