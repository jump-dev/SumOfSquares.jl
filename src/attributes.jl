"""
    struct CertificateBasis <: MOI.AbstractConstraintAttribute end

A constraint attribute for the basis indexing the
[`GramMatrixAttribute`](@ref) and [`MomentMatrixAttribute`](@ref) certificates.
"""
struct CertificateBasis <: MOI.AbstractConstraintAttribute end

"""
    GramMatrixAttribute(
        multiplier_index::Int = 0,
        result_index::Int = 1,
    )

A constraint attribute for the [`GramMatrix`](@ref) of the `multiplier_index`th
Sum-of-Squares polynomial ``s_i(x)`` where ``i`` is `multiplier_index` of the certificate:
```math
p(x) = s_0(x) + w_1(x) s_1(x) + \\cdots + w_m(x) s_m(x)
```
The gram matrix of a Sum-of-Squares polynomial ``s_i(x)`` is the
the positive semidefinite matrix ``Q`` such that ``s_i(x) = b_i(x)^\\top Q b_i(x)``
where ``b_i(x)`` is the gram basis.
"""
@kwdef struct GramMatrixAttribute <: MOI.AbstractConstraintAttribute
    multiplier_index::Int = 0
    result_index::Int = 1
end

"""
    @kwdef struct SOSDecompositionAttribute
        multiplier_index::Int = 0
        ranktol::Real
        dec::MultivariateMoments.LowRankLDLTAlgorithm
        result_index::Int
    end

A constraint attribute for the [`SOSDecomposition`](@ref) of a constraint.
By default, it is computed using
`SOSDecomposition(gram, ranktol, dec)` where `gram` is the value of the
[`GramMatrixAttribute`](@ref).
"""
@kwdef struct SOSDecompositionAttribute <: MOI.AbstractConstraintAttribute
    multiplier_index::Int = 0
    ranktol::Real
    dec::MultivariateMoments.LowRankLDLTAlgorithm
    result_index::Int = 1
end

function MOI.get_fallback(
    model::MOI.ModelLike,
    attr::SOSDecompositionAttribute,
    ci::MOI.ConstraintIndex,
)
    gram = MOI.get(model, attr.result_index, ci)
    return SOSDecomposition(gram, attr.ranktol, attr.dec)
end

"""
    MomentMatrixAttribute(
        multiplier_index::Int = 0,
        result_index::Int = 1,
    )

A constraint attribute for the `MomentMatrix` of the `multiplier_index`th
Sum-of-Squares polynomial ``s_i(x)`` where ``i`` is `multiplier_index` of the certificate:
```math
p(x) = s_0(x) + w_1(x) s_1(x) + \\cdots + w_m(x) s_m(x)
```
It corresponds to the dual of the Sum-of-Squares constraint for the constraint
for ``s_i(x)`` to be a Sum-of-Squares.
"""
@kwdef struct MomentMatrixAttribute <: MOI.AbstractConstraintAttribute
    multiplier_index::Int = 0
    result_index::Int = 1
end

"""
    struct MultiplierIndexBoundsError{AttrType} <: Exception
        attr::AttrType
        range::UnitRange{Int}
    end

An error indicating that the requested attribute `attr` could not be retrieved,
because the multiplier index is out of the range of valid indices.
"""
struct MultiplierIndexBoundsError{AttrType} <: Exception
    attr::AttrType
    range::UnitRange{Int}
end

function check_multiplier_index_bounds(attr, range)
    if !(attr.multiplier_index in range)
        throw(MultiplierIndexBoundsError(attr, range))
    end
end

function Base.showerror(io::IO, err::MultiplierIndexBoundsError)
    return print(
        io,
        "Multiplier index of attribute $(err.attr) out of bounds. The index " *
        "must be in the range $(err.range).",
    )
end

"""
    LagrangianMultipliers(result_index::Int)
    LagrangianMultipliers()

A constraint attribute fot the `LagrangianMultipliers` associated to the
inequalities of the domain of a constraint. There is one multiplier per
inequality and no multiplier for equalities as the equalities are handled by
reducing the polynomials over the ideal they generate instead of explicitely
creating multipliers.
"""
struct LagrangianMultipliers <: MOI.AbstractConstraintAttribute
    result_index::Int
end
LagrangianMultipliers() = LagrangianMultipliers(1)

const _Attributes = Union{
    CertificateBasis,
    GramMatrixAttribute,
    SOSDecompositionAttribute,
    MomentMatrixAttribute,
    LagrangianMultipliers,
}

# Needs to declare it set by optimize that it is not queried in the Caching
# optimize, even of `CertificateBasis` which is set before optimize.
MOI.is_set_by_optimize(::_Attributes) = true

# If a variable is bridged, the `VectorOfVariables`-in-`SOSPolynomialSet` is
# bridged by `MOI.Bridges.Constraint.VectorFunctionizeBridge` and it has
# to pass the constraint to the SOS bridge.
MOI.Bridges.Constraint.invariant_under_function_conversion(::_Attributes) = true

# They do not contain any variable so `substitute_variables` would be the identity.
# We cannot just implement `substitute_variables` since it some variables cannot
# be unbridged.
function MOI.Bridges.unbridged_function(
    ::MOI.Bridges.AbstractBridgeOptimizer,
    value::Union{GramMatrix{T},MultivariateMoments.MomentVector{T}},
) where {T<:Number}
    return value
end

# This is type piracy but we tolerate it.
const ObjectWithoutIndex = Union{
    AbstractGramMatrix{<:MOI.Utilities.ObjectWithoutIndex},
    SOSDecomposition{<:MOI.Utilities.ObjectWithoutIndex},
}
const ObjectOrTupleWithoutIndex =
    Union{ObjectWithoutIndex,Tuple{Vararg{ObjectWithoutIndex}}}
const ObjectOrTupleOrArrayWithoutIndex =
    Union{ObjectOrTupleWithoutIndex,AbstractArray{<:ObjectOrTupleWithoutIndex}}
MOI.Utilities.map_indices(::Function, x::ObjectOrTupleOrArrayWithoutIndex) = x
function MOI.Utilities.substitute_variables(
    ::Function,
    x::ObjectOrTupleOrArrayWithoutIndex,
)
    return x
end
