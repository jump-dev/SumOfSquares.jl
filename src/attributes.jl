struct CertificateMonomials <: MOI.AbstractConstraintAttribute end
struct GramMatrixAttribute <: MOI.AbstractConstraintAttribute end
struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute end
struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute end

function MOI.is_set_by_optimize(::Union{CertificateMonomials,
                                        GramMatrixAttribute,
                                        MomentMatrixAttribute,
                                        LagrangianMultipliers})
    return true
end
