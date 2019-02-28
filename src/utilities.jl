# + is not defined between SingleVariable
function MP.polynomial(p::GramMatrix{MOI.SingleVariable})
    Q = convert(Vector{MOI.ScalarAffineFunction{Float64}}, p.Q.Q)
    MP.polynomial(GramMatrix(SymMatrix(Q, p.Q.n), p.x))
end
function MP.polynomial(p::GramMatrix{F}) where {F <: MOI.AbstractFunction}
    MP.polynomial(p, MOIU.promote_operation(+, Float64, F, F))
end

function primal_value(model, p::GramMatrix{MOI.SingleVariable})
    # TODO [perf] use MOI typed mapped array
    Q = MOI.get(model, MOI.VariablePrimal(),
                MOI.VariableIndex[sv.variable for sv in p.Q.Q])
    return GramMatrix(SymMatrix(Q, p.Q.n), p.x)
end

function add_matrix_variable_bridge(
    model::MOI.ModelLike, MCT,
    side_dimension::Integer, T::Type)
    mat_cone = matrix_cone(MCT, side_dimension)
    VB = variable_bridge_type(typeof(mat_cone), T)
    return add_variable_bridge(VB, model, mat_cone)
end
function union_vector_bridge_types(
    MCT, T::Type)
    return Union{variable_bridge_type(typeof(matrix_cone(MCT, 0)), T),
                 variable_bridge_type(typeof(matrix_cone(MCT, 1)), T),
                 variable_bridge_type(typeof(matrix_cone(MCT, 2)), T),
                 variable_bridge_type(typeof(matrix_cone(MCT, 3)), T)}
end
function append_added_constraint_types(
    added, MCT, T::Type)
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 0)), T)))
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 1)), T)))
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 2)), T)))
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 3)), T)))
    return added
end
