# `+` is not defined between `MOI.SingleVariable`.
_promote_sum(::Type{MOI.SingleVariable}) = MOI.ScalarAffineFunction{Float64}
_promote_add_mul(::Type{MOI.SingleVariable}) = MOI.ScalarQuadraticFunction{Float64}
function MP.polynomial(p::GramMatrix{MOI.SingleVariable})
    Q = convert(Vector{MOI.ScalarAffineFunction{Float64}}, p.Q.Q)
    return MP.polynomial(GramMatrix(SymMatrix(Q, p.Q.n), p.x))
end
#function MP.polynomial(p::GramMatrix{F}) where {F <: MOI.AbstractFunction}
#    MP.polynomial(p, MOIU.promote_operation(+, Float64, F, F))
#end

function primal_value(model, p::GramMatrix{MOI.SingleVariable})
    # TODO [perf] use MOI typed mapped array
    Q = MOI.get(model, MOI.VariablePrimal(),
                MOI.VariableIndex[sv.variable for sv in p.Q.Q])
    return GramMatrix(SymMatrix(Q, p.Q.n), p.x)
end
