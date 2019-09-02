function union_vector_types(MCT, T::Type)
    return Union{MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 0))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 1))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 2))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 3))}}
end
