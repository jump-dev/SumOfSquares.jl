struct CertificateMonomials <: MOI.AbstractConstraintAttribute end
struct GramMatrix <: MOI.AbstractConstraintAttribute end
struct MomentMatrix <: MOI.AbstractConstraintAttribute end
struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute end

function MOI.is_set_by_optimize(::Union{CertificateMonomials, GramMatrix,
                                        MomentMatrix, LagrangianMultipliers})
    return true
end
