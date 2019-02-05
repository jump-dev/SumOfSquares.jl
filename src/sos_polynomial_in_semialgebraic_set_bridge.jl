function lagrangian_multiplier(model::MOI.ModelLike, p, set::SOSSubCones, q, mindegree::Integer, maxdegree::Integer)
    mindegree_q, maxdegree_q = extdegree(q)
    # extdegree's that s^2 should have so that s^2 * p has degrees between mindegree and maxdegree
    mindegree_s2 = mindegree - mindegree_q
    maxdegree_s2 = maxdegree - maxdegree_q
    # extdegree's for s
    mindegree_s = max(0, div(mindegree_s2, 2))
    # If maxdegree_s2 is odd, div(maxdegree_s2,2) would make s^2 have degree up to maxdegree_s2-1
    # for this reason, we take div(maxdegree_s2+1,2) so that s^2 have degree up to maxdegree_s2+1
    maxdegree_s = div(maxdegree_s2 + 1, 2)
    # FIXME handle the case where `p`, `q_i`, ...  do not have the same variables
    # so instead of `variable(p)` we would have the union of them all
    @assert variables(q) ⊆ variables(p)
    monos = monomials(variables(p), mindegree_s:maxdegree_s)
    return matpoly_in_cone(model, monos, set)
end

struct SOSPolynomialInSemialgebraicSetBridge{T, F <: MOI.AbstractVectorFunction,
                                             DT <: AbstractSemialgebraicSet,
                                             CT <: SOSSubCones,
                                             BT <: PolyJuMP.AbstractPolynomialBasis,
                                             MT <: AbstractMonomial,
                                             MVT <: AbstractVector{MT},
                                             NPT <: Tuple} <: MOIB.AbstractBridge
    lagrangian_multipliers::Vector{MatPolynomial{MOI.SingleVariable, MT, MVT}}
    constraint::MOI.ConstraintIndex{F, SOSPolynomialSet{DT, CT, BT, MT, MVT, NPT}}
end

function SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, BT, MT, MVT, NPT}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    set::SOSPolynomialSet{<:BasicSemialgebraicSet}) where {T, F, DT, CT, BT, MT, MVT, NPT}

    @assert MOI.output_dimension(f) == length(set.monomials)
    p = polynomial(collect(MOIU.eachscalar(f)), set.monomials)
    λ = lagrangian_multiplier.(model, p, set.cone, set.domain.p, set.mindegree,
                               set.maxdegree)
    p -= dot(set.domain.p, λ)
    monos = monomials(p)
    new_set = SOSPolynomialSet(set.domain.V, set.cone, set.basis, monos,
                               set.newton_polytope, set.mindegree,
                               set.maxdegree)
    constraint = MOI.add_constraint(model, MOIU.vectorize(coefficients(p)),
                                    new_set)

    return SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, BT, MT, MVT, NPT}(λ, constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialInSemialgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{<:BasicSemialgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, BT, MT, MVT, NPT}}) where {T, F, DT, CT, BT, MT, MVT, NPT}
    return [(F, SOSPolynomialSet{DT, CT, BT, MT, MVT, NPT})]
end
function MOIB.concrete_bridge_type(::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:SOSPolynomialSet{BasicSemialgebraicSet{S, PS, AT}, CT, BT, MT, MVT, NPT}}) where {T, S, PS, AT, CT, BT, MT, MVT, NPT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    return SOSPolynomialInSemialgebraicSetBridge{T, G, AT, CT, BT, MT, MVT, NPT}
end

# Attributes, Bridge acting as an model
function MOI.get(::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT, NPT},
                 ::MOI.NumberOfConstraints{F, SOSPolynomialSet{DT, BT, MT, MVT, NPT}}) where {T, F, DT, BT, MT, MVT, NPT}
    return 1
end
function MOI.get(b::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT, NPT},
                 ::MOI.ListOfConstraintIndices{F, SOSPolynomialSet{DT, BT, MT, MVT, NPT}}) where {T, F, DT, BT, MT, MVT, NPT}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::SOSPolynomialInSemialgebraicSetBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint

# The monomials might be different from the ones of the original polynomial
# because of the ∑ λ_i s_i(x) so we don't define ConstraintPrimal and
# ConstraintDual, as the caller won't know how to reshape it
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MomentMatrix, CertificateMonomials},
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, bridge.constraint)
end
function MOI.get(model::MOI.ModelLike, ::LagrangianMultipliers,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return map(λ -> primal_value(model, λ), bridge.lagrangian_multipliers)
end
