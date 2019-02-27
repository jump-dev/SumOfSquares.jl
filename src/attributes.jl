"""
    struct CertificateMonomials <: MOI.AbstractConstraintAttribute end

A constraint attribute for the monomials indexing the
[`GramMatrixAttribute`](@ref) and [`MomentMatrixAttribute`](@ref) certificates.
"""
struct CertificateMonomials <: MOI.AbstractConstraintAttribute end

"""
    struct GramMatrixAttribute <: MOI.AbstractConstraintAttribute end

A constraint attribute for the [`GramMatrix`](@ref) of a constraint, that is,
the positive semidefinte matrix `Q` indexed by the monomials in the vector `X`
such that ``X^\\top Q X`` is the sum-of-squares certificate of the constraint.
The
"""
struct GramMatrixAttribute <: MOI.AbstractConstraintAttribute end

"""
    struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute end

A constraint attribute fot the `MomentMatrix` of a constraint.
"""
struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute end

"""
    struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute end

A constraint attribute fot the `LagrangianMultipliers` assiciated to the
inequalities of the domain of a constraint. There is one multiplier per
inequality and no multiplier for equalities as the equalities are handled by
reducing the polynomials over the ideal they generate instead of explicitely
creating multipliers.
"""
struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute end

function MOI.is_set_by_optimize(::Union{CertificateMonomials,
                                        GramMatrixAttribute,
                                        MomentMatrixAttribute,
                                        LagrangianMultipliers})
    return true
end
