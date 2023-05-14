function union_constraint_types(MCT)
    UMCT = SOS.union_constraint_indices_types(MCT)
    return Union{UMCT,Vector{UMCT}}
    # `Vector{UMCT}` is for sparse SOS
    # As basis can have different sizes, we need `Vector{Union{...}}` instead of
    # `Union{Vector, ...}`.
end
function union_set_types(MCT)
    return Union{
        SOS.matrix_cone_type(MCT, 0),
        SOS.matrix_cone_type(MCT, 1),
        SOS.matrix_cone_type(MCT, 2),
        SOS.matrix_cone_type(MCT, 3),
    }
end
function constrained_variable_types(MCT)
    return [
        (SOS.matrix_cone_type(MCT, 0),),
        (SOS.matrix_cone_type(MCT, 1),),
        (SOS.matrix_cone_type(MCT, 2),),
        (SOS.matrix_cone_type(MCT, 3),),
    ]
end
