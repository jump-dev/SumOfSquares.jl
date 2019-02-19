function lagrangian_multiplier(model::MOI.ModelLike, p, set::SOSLikeCone, q,
                               mindegree::Integer, maxdegree::Integer, VB::Type)
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
    Q, variable_bridge = add_variable_bridge(
        VB, model, matrix_cone(matrix_cone_type(typeof(set)), length(monos)))
    return build_gram_matrix(Q, monos), variable_bridge, monos
end

struct SOSPolynomialInSemialgebraicSetBridge{
    T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
    CT <: SOSLikeCone, VB <: AbstractVariableBridge,
    BT <: PolyJuMP.AbstractPolynomialBasis, MT <: AbstractMonomial,
    MVT <: AbstractVector{MT}, NPT <: Tuple} <: MOIB.AbstractBridge
    lagrangian_monomials::Vector{MVT}
    lagrangian_bridges::Vector{VB}
    constraint::MOI.ConstraintIndex{F, SOSPolynomialSet{DT, CT, BT, MT, MVT, NPT}}
end

function SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, VB, BT, MT, MVT, NPT}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    set::SOSPolynomialSet{<:BasicSemialgebraicSet}) where {
        # Need to specify types to avoid ambiguity with the default constructor
        T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
        CT <: SOSLikeCone, VB <: AbstractVariableBridge,
        BT <: PolyJuMP.AbstractPolynomialBasis, MT <: AbstractMonomial,
        MVT <: AbstractVector{MT}, NPT <: Tuple
    }
    @assert MOI.output_dimension(f) == length(set.monomials)
    p = polynomial(collect(MOIU.eachscalar(f)), set.monomials)
    n = length(set.domain.p)
    λ_monos   = Vector{MVT}(undef, n)
    λ_bridges = Vector{VB}(undef, n)
    for (i, q) in enumerate(set.domain.p)
        λ, λ_bridges[i], λ_monos[i] = lagrangian_multiplier(
            model, p, set.cone, q, set.mindegree, set.maxdegree, VB)
        p -= λ * q
    end
    monos = monomials(p)
    new_set = SOSPolynomialSet(set.domain.V, set.cone, set.basis, monos,
                               set.newton_polytope, set.mindegree,
                               set.maxdegree)
    constraint = MOI.add_constraint(model, MOIU.vectorize(coefficients(p)),
                                    new_set)

    return SOSPolynomialInSemialgebraicSetBridge{
        T, F, DT, CT, VB, BT, MT, MVT, NPT}(λ_monos, λ_bridges, constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialInSemialgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{<:BasicSemialgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, VB, BT, MT, MVT, NPT}}) where {T, F, DT, CT, VB, BT, MT, MVT, NPT}
    added = [(F, SOSPolynomialSet{DT, CT, BT, MT, MVT, NPT})]
    append!(added, MOIB.added_constraint_types(VB))
    return added
end
function MOIB.concrete_bridge_type(::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:SOSPolynomialSet{BasicSemialgebraicSet{S, PS, AT}, CT, BT, MT, MVT, NPT}}) where {T, S, PS, AT, CT, BT, MT, MVT, NPT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    VB = variable_bridge_type(matrix_cone_type(CT), T)
    return SOSPolynomialInSemialgebraicSetBridge{T, G, AT, CT, VB, BT, MT, MVT, NPT}
end

# Attributes, Bridge acting as an model
function MOI.get(::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, VB, BT, MT, MVT, NPT},
                 ::MOI.NumberOfConstraints{F, SOSPolynomialSet{DT, BT, MT, MVT, NPT}}) where {T, F, DT, VB, BT, MT, MVT, NPT}
    return 1
end
function MOI.get(b::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, VB, BT, MT, MVT, NPT},
                 ::MOI.ListOfConstraintIndices{F, SOSPolynomialSet{DT, BT, MT, MVT, NPT}}) where {T, F, DT, VB, BT, MT, MVT, NPT}
    return [b.constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::SOSPolynomialInSemialgebraicSetBridge)
    MOI.delete(model, bridge.constraint)
    for variable_bridge in bridge.lagrangian_bridges
        MOI.delete(model, variable_bridge)
    end
end

# Attributes, Bridge acting as a constraint

# The monomials might be different from the ones of the original polynomial
# because of the ∑ λ_i s_i(x) so we don't define ConstraintPrimal and
# ConstraintDual, as the caller won't know how to reshape it
function MOI.get(model::MOI.ModelLike,
                 attr::Union{CertificateMonomials, GramMatrixAttribute,
                             MomentMatrixAttribute},
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, bridge.constraint)
end
function MOI.get(model::MOI.ModelLike, ::LagrangianMultipliers,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    @assert eachindex(bridge.lagrangian_bridges) == eachindex(bridge.lagrangian_monomials)
    map(i -> build_gram_matrix(MOI.get(model, GramMatrixAttribute(),
                                       bridge.lagrangian_bridges[i]),
                               bridge.lagrangian_monomials[i]),
        eachindex(bridge.lagrangian_bridges))
end
