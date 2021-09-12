"""
    struct CertificateBasis <: MOI.AbstractConstraintAttribute end

A constraint attribute for the basis indexing the
[`GramMatrixAttribute`](@ref) and [`MomentMatrixAttribute`](@ref) certificates.
"""
struct CertificateBasis <: MOI.AbstractConstraintAttribute end

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
    struct SOSDecompositionAttribute
        ranktol::Real
        dec::MultivariateMoments.LowRankChol
    end

A constraint attribute for the [`SOSDecomposition`](@ref) of a constraint.
By default, it is computed using
`SOSDecomposition(gram, ranktol, dec)` where `gram` is the value of the
[`GramMatrixAttribute`](@ref).
"""
struct SOSDecompositionAttribute <: MOI.AbstractConstraintAttribute
    ranktol::Real
    dec::MultivariateMoments.LowRankChol
    result_index::Int
end
function SOSDecompositionAttribute(ranktol::Real, dec::MultivariateMoments.LowRankChol)
    return SOSDecompositionAttribute(ranktol, dec, 1)
end

function MOI.get_fallback(
    model::MOI.ModelLike,
    attr::SOSDecompositionAttribute,
    ci::MOI.ConstraintIndex,
)
    gram = MOI.get(model, attr.result_index, ci)
    return SOSDecomposition(gram, attr.ranktol, attr.dec)
end
# TODO bridges should redirect to `MOI.get_fallback` as well so that
# we can just use `Union{MOI.ConstraintIndex,MOI.Bridges.AbstractBridge}` above:
function MOI.get(
    model::MOI.ModelLike,
    attr::SOSDecompositionAttribute,
    bridge::Union{Bridges.Constraint.SOSPolynomialBridge,Bridges.Constraint.SOSPolynomialInSemialgebraicSetBridge},
)
    gram = MOI.get(model, GramMatrixAttribute(attr.result_index), bridge)
    return SOSDecomposition(gram, attr.ranktol, attr.dec)
end

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
# optimize, even of `CertificateBasis` which is set befor optimize.
function MOI.is_set_by_optimize(::Union{CertificateBasis,
                                        GramMatrixAttribute,
                                        SOSDecompositionAttribute,
                                        MomentMatrixAttribute,
                                        LagrangianMultipliers})
    return true
end

# If a variable is bridged, the `VectorOfVariables`-in-`SOSPolynomialSet` is
# bridged by `MOI.Bridges.Constraint.VectorFunctionizeBridge` and it has
# to pass the constraint to the SOS bridge.
function MOI.Bridges.Constraint.invariant_under_function_conversion(::Union{
    CertificateBasis, GramMatrixAttribute,
    MomentMatrixAttribute, LagrangianMultipliers})
    return true
end

# This is type piracy but we tolerate it.
const ObjectWithoutIndex = Union{AbstractGramMatrix{<:MOI.Utilities.ObjectWithoutIndex},SOSDecomposition{<:MOI.Utilities.ObjectWithoutIndex}}
const ObjectOrTupleWithoutIndex = Union{ObjectWithoutIndex, Tuple{Vararg{ObjectWithoutIndex}}}
const ObjectOrTupleOrArrayWithoutIndex = Union{ObjectOrTupleWithoutIndex, AbstractArray{<:ObjectOrTupleWithoutIndex}}
MOI.Utilities.map_indices(::Function, x::ObjectOrTupleOrArrayWithoutIndex) = x
MOI.Utilities.substitute_variables(::Function, x::ObjectOrTupleOrArrayWithoutIndex) = x
