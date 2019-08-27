"""
    struct CertificateMonomials <: MOI.AbstractConstraintAttribute end

A constraint attribute for the monomials indexing the
[`GramMatrixAttribute`](@ref) and [`MomentMatrixAttribute`](@ref) certificates.
"""
struct CertificateMonomials <: MOI.AbstractConstraintAttribute end

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

"""
    MomentMatrixAttribute(N)
    MomentMatrixAttribute()

A constraint attribute fot the `MomentMatrix` of a constraint.
"""
struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute
    N::Int
end
MomentMatrixAttribute() = MomentMatrixAttribute(1)

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

# This is type piracy but we tolerate it.
const ObjectWithoutIndex = Union{MultivariateMoments.MomentMatrix{<:MOI.Utilities.ObjectWithoutIndex}, GramMatrix{<:MOI.Utilities.ObjectWithoutIndex}}
const ObjectOrTupleWithoutIndex = Union{ObjectWithoutIndex, Tuple{Vararg{ObjectWithoutIndex}}}
const ObjectOrTupleOrArrayWithoutIndex = Union{ObjectOrTupleWithoutIndex, AbstractArray{<:ObjectOrTupleWithoutIndex}}
MOI.Utilities.map_indices(::Function, x::ObjectOrTupleOrArrayWithoutIndex) = x
MOI.Utilities.substitute_variables(::Function, x::ObjectOrTupleOrArrayWithoutIndex) = x
