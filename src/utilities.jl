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

### Utilities for writting code that works but at the JuMP and MOI level ###

function _add_constraint(model::MOI.ModelLike, f::Vector{MOI.SingleVariable}, s)
    return _add_constraint(model, MOI.VariableIndex[v.variable for v in f], s)
end
#function _add_constraint(model::MOI.ModelLike, f::Vector{<:MOI.AbstractScalarFunction}, s)
#    return _add_constraint(model, MOIU.vectorize(f), s)
#end
function _add_constraint(model::MOI.ModelLike, f, s)
    MOI.add_constraint(model, f, s)
end
function _add_constraint(model::JuMP.AbstractModel, f, s)
    JuMP.add_constraint(model, JuMP.build_constraint(error, f, s))
end

_vcat(funs::MOI.AbstractFunction...) = MOIU.operate(vcat, Float64, funs...)
_vcat(funs::JuMP.AbstractJuMPScalar...) = vcat(funs...)

### Utilities for constraining the symmetric matrix of a MatPolynomial ###

# PSD constraints on 2x2 matrices are SOC representable.
# [Q11 Q12] is PSD iff Q11, Q22 ≥ 0 and       Q11*Q22 ≥     Q12 ^2
# [Q12 Q22] is PSD                      <=> 2*Q11*Q22 ≥ (√2*Q12)^2
function soc_psd_constraint(model, Q11, Q12, Q22)
    _add_constraint(model, _vcat(Q11, Q22, √2*Q12),
                    MOI.RotatedSecondOrderCone(3))
end

function matrix_add_constraint(model, Q::Vector, set::MOI.AbstractVectorSet)
    @assert length(Q) == MOI.dimension(set)
    n = set.side_dimension
    if n == 1
        # PSD constraints on 1x1 matrices are equivalent to the
        # nonnegativity of the only entry
        return _add_constraint(model, Q[1], MOI.GreaterThan(0.0))
    elseif false && n == 2 && !(set isa DiagonallyDominantConeTriangle)
        # PSD constraints on 2x2 matrices are SOC representable.
        # For the DD cone, we want to avoid using SOC and only use LP
        # TODO check that ConstraintPrimal and ConstraintDual are not used
        #      we should probably make a bridge that transform them
        return soc_psd_constraint(model, Q...)
        # FIXME for SDOI solvers, it will be transformed to a 3x3 matrix as they
        #       do not support SOC.
        #       We should transform it to RSOC and make an particular case for RSOC->PSD bridge of dimension 3. This exception would also help for the geometric cone which also generates dimension 3 RSOC
    else
        # PSD constraints on nxn matrices with n ≥ 3 is not SOC representable,
        # see [F18].
        #
        # [F18] Fawzi, Hamza
        # On representing the positive semidefinite cone using the second-order
        # cone.
        # Mathematical Programming (2018): 1-10.
        return _add_constraint(model, Q, set)
    end
end
function matrix_add_constraint(model, p::MatPolynomial, S::Type)
    matrix_add_constraint(model, p.Q.Q, S(length(p.x)))
end
