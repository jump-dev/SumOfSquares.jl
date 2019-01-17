### Utilities for writting code that works but at the JuMP and MOI level ###

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
function soc_psd_constraint(model, Q11, Q12, Q22)
    _add_constraint(model, _vcat(Q11 + Q22, 2.0Q12, Q11 - Q22),
                    MOI.SecondOrderCone(3))
end

function matrix_add_constraint(model, Q::Vector, set::MOI.AbstractVectorSet)
    @assert length(Q) == MOI.dimension(set)
    n = set.side_dimension
    if n == 1
        # PSD constraints on 1x1 matrices are equivalent to the
        # nonnegativity of the only entry
        _add_constraint(model, Q[1], MOI.GreaterThan(0.0))
    elseif n == 2 && !(S <: DiagonallyDominantConeTriangle)
        # PSD constraints on 2x2 matrices are SOC representable.
        # For the DD cone, we want to avoid using SOC and only use LP
        soc_psd_constraint(model, Q...)
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
        _add_constraint(model, Q, set)
    end
end
function matrix_add_constraint(model, p::MatPolynomial, S::Type)
    matrix_add_constraint(model, p.Q.Q, S(length(p.x)))
end
