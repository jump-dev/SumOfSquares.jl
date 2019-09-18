function union_constraint_types(MCT)
    return Union{MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 0))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 1))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 2))},
                 MOI.ConstraintIndex{MOI.VectorOfVariables, typeof(SOS.matrix_cone(MCT, 3))}}
end
function union_set_types(MCT)
    return Union{typeof(SOS.matrix_cone(MCT, 0)),
                 typeof(SOS.matrix_cone(MCT, 1)),
                 typeof(SOS.matrix_cone(MCT, 2)),
                 typeof(SOS.matrix_cone(MCT, 3))}
end
function constrained_variable_types(MCT)
    return [
        (typeof(SOS.matrix_cone(MCT, 0)),),
        (typeof(SOS.matrix_cone(MCT, 1)),),
        (typeof(SOS.matrix_cone(MCT, 2)),),
        (typeof(SOS.matrix_cone(MCT, 3)),)
   ]
end
