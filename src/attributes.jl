"""
    struct CertificateMonomials <: MOI.AbstractConstraintAttribute end

A constraint attribute for the monomials indexing the
[`GramMatrixAttribute`](@ref) and [`MomentMatrixAttribute`](@ref) certificates.
"""
struct CertificateMonomials <: MOI.AbstractConstraintAttribute end
# This is type piracy but we tolerate it.
MOIU.map_indices(::Function, monovec::AbstractVector{<:MP.AbstractMonomial}) = monovec
MOIU.substitute_variables(::Function, monovec::AbstractVector{<:MP.AbstractMonomial}) = monovec

"""
    GramMatrixAttribute(N)
    GramMatrixAttribute()

A constraint attribute for the [`GramMatrix`](@ref) of a constraint, that is,
the positive semidefinite matrix `Q` indexed by the monomials in the vector `X`
such that ``X^\\top Q X`` is the sum-of-squares certificate of the constraint.
"""
struct GramMatrixAttribute <: MOI.AbstractConstraintAttribute
    N::Int
end
GramMatrixAttribute() = GramMatrixAttribute(1)
MOIU.map_indices(::Function, gram::GramMatrix{<:MOIU.ObjectWithoutIndex}) = gram
MOIU.substitute_variables(::Function, gram::GramMatrix{<:MOIU.ObjectWithoutIndex}) = gram

"""
    MomentMatrixAttribute(N)
    MomentMatrixAttribute()

A constraint attribute fot the `MomentMatrix` of a constraint.
"""
struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute
    N::Int
end
MomentMatrixAttribute() = MomentMatrixAttribute(1)
# This is type piracy but we tolerate it.
MOIU.map_indices(::Function, mom::MultivariateMoments.MomentMatrix{<:MOIU.ObjectWithoutIndex}) = mom
MOIU.substitute_variables(::Function, mom::MultivariateMoments.MomentMatrix{<:MOIU.ObjectWithoutIndex}) = mom

"""
    LagrangianMultipliers(N)
    LagrangianMultipliers()

A constraint attribute fot the `LagrangianMultipliers` associated to the
inequalities of the domain of a constraint. There is one multiplier per
inequality and no multiplier for equalities as the equalities are handled by
reducing the polynomials over the ideal they generate instead of explicitely
creating multipliers.
"""
struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute
    N::Int
end
LagrangianMultipliers() = LagrangianMultipliers(1)

# Needs to declare it set by optimize that it is not queried in the Caching
# optimize, even of `CertificateMonomials` which is set befor optimize.
function MOI.is_set_by_optimize(::Union{CertificateMonomials,
                                        GramMatrixAttribute,
                                        MomentMatrixAttribute,
                                        LagrangianMultipliers})
    return true
end
