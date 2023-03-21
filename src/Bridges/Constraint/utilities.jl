function union_constraint_types(MCT)
    UMCT = SOS.union_constraint_indices_types(MCT)
    return Union{UMCT,Vector{UMCT}}
    # `Vector{UMCT}` is for sparse SOS
    # As basis can have different sizes, we need `Vector{Union{...}}` instead of
    # `Union{Vector, ...}`.
end
function union_set_types(MCT)
    return Union{
        typeof(SOS.matrix_cone(MCT, 0)),
        typeof(SOS.matrix_cone(MCT, 1)),
        typeof(SOS.matrix_cone(MCT, 2)),
        typeof(SOS.matrix_cone(MCT, 3)),
    }
end
function constrained_variable_types(MCT)
    return [
        (typeof(SOS.matrix_cone(MCT, 0)),),
        (typeof(SOS.matrix_cone(MCT, 1)),),
        (typeof(SOS.matrix_cone(MCT, 2)),),
        (typeof(SOS.matrix_cone(MCT, 3)),),
    ]
end
