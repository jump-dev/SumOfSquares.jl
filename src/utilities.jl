# + is not defined between SingleVariable
function MP.polynomial(p::MatPolynomial{MOI.SingleVariable})
    Q = convert(Vector{MOI.ScalarAffineFunction{Float64}}, p.Q.Q)
    MP.polynomial(MatPolynomial(SymMatrix(Q, p.Q.n), p.x))
end
function MP.polynomial(p::MatPolynomial{F}) where {F <: MOI.AbstractFunction}
    MP.polynomial(p, MOIU.promote_operation(+, Float64, F, F))
end

function primal_value(model, p::MatPolynomial{MOI.SingleVariable})
    # TODO [perf] use MOI typed mapped array
    Q = MOI.get(model, MOI.VariablePrimal(),
                MOI.VariableIndex[sv.variable for sv in p.Q.Q])
    return MatPolynomial(SymMatrix(Q, p.Q.n), p.x)
end

function add_matrix_variable_bridge(
    model::MOI.ModelLike, MCT::Type{<:MOI.AbstractVectorSet},
    side_dimension::Integer, T::Type)
    mat_cone = matrix_cone(MCT, side_dimension)
    VB = variable_bridge_type(typeof(mat_cone), T)
    return add_variable_bridge(VB, model, mat_cone)
end
function union_vector_bridge_types(
    MCT::Type{<:MOI.AbstractVectorSet}, T::Type)
    return Union{variable_bridge_type(typeof(matrix_cone(MCT, 1)), T),
                 variable_bridge_type(typeof(matrix_cone(MCT, 2)), T),
                 variable_bridge_type(typeof(matrix_cone(MCT, 3)), T)}
end
function append_added_constraint_types(
    added, MCT::Type{<:MOI.AbstractVectorSet}, T::Type)
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 1)), T)))
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 2)), T)))
    append!(added, MOIB.added_constraint_types(variable_bridge_type(typeof(matrix_cone(MCT, 3)), T)))
    return added
end
