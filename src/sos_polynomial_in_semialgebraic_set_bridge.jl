function lagrangian_multiplier(model::JuMP.Model, p, set::SOSSubCones, q, mindegree::Integer, maxdegree::Integer)
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
    return JuMP.add_variable(model, PolyJuMP.Variable(_varconetype(set)(monos),
                                                      false, false))
end

struct SOSPolynomialInSemialgebraicSetBridge{T, F <: MOI.AbstractVectorFunction,
                                             DT <: AbstractSemialgebraicSet,
                                             BT <: PolyJuMP.AbstractPolynomialBasis,
                                             MT <: AbstractMonomial,
                                             MVT <: AbstractVector{MT},
                                             NPT <: Tuple} <: MOIB.AbstractBridge
    constraint::MOI.ConstraintIndex{F, MOI.SOSPolynomialSet{DT, BT, MT, MVT, NPT}}
end

function SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT, NPT}(model::MOI.ModelLike,
                                                                           f::MOI.AbstractVectorFunction,
                                                                           s::SOSPolynomialInSemialgebraicSetBridge{<:AbstractSemialgebraicSet}) where {T, F, DT, BT, MT, MVT, NPT}
    @assert MOI.output_dimension(f) == length(s.monomials)
    p = polynomial(MOIU.eachscalar(f), s.monomials)
    λ = lagrangian_multiplier.(model, p, Ref(set), domain.p, mindegree, maxdegree)
    p -= dot(λ, domain.p)
    constraint = PolyJuMP.addpolyconstraint!(m, p, set, domain.V, basis; kws...)
    constraint.lagrangian_multipliers = λ

    return SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT}(zero_constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialInSemialgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{FullSpace}}) where T
    return true
end
function added_constraint_types(::Type{SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT}}) where {T, F}
    return [(F, MOI.ZeroPolynomialSet{DT, BT, MT, MVT})]
end
function concrete_bridge_type(::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              ::Type{SOSPolynomialSet{FullSpace, <:MonomialBasis, MT, MVT}}) where {T, MT, MVT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, zeros(T))
    return SOSPolynomialInSemialgebraicSetBridge{T, G, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, MOI.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return 1
end
function MOI.get(b::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, BT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, MOI.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::SOSPolynomialInSemialgebraicSetBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 c::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, c.zero_constraint)
end
