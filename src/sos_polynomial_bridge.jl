struct SOSPolynomialBridge{T, F <: MOI.AbstractVectorFunction,
                           DT <: AbstractSemialgebraicSet,
                           BT <: PolyJuMP.AbstractPolynomialBasis,
                           MT <: AbstractMonomial,
                           MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    gram_matrix::MatPolynomial{MOI.SingleVariable, MT, MVT}
    gram_constraint::MOI.ConstraintIndex{MOI.VectorOfVariables} # TODO add set type
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}
end

function SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::SOSPolynomialSet{<:AbstractAlgebraicSet}) where {
        # Need to specify types to avoid ambiguity with the default constructor
        T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
        BT <: PolyJuMP.AbstractPolynomialBasis, MT <: AbstractMonomial,
        MVT <: AbstractVector{MT}
    }
    @assert MOI.output_dimension(f) == length(s.monomials)
    p = polynomial(collect(MOIU.eachscalar(f)), s.monomials)
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = convert(typeof(p), rem(p, ideal(s.domain)))
    X = monomials_half_newton_polytope(monomials(r), s.newton_polytope)
    g, gram_matrix, gram_constraint = gram_in_cone(model, X, s.cone, T)
    q = r - g
    set = PolyJuMP.ZeroPolynomialSet(s.domain, s.basis, monomials(q))
    coefs = MOIU.vectorize(coefficients(q))
    zero_constraint = MOI.add_constraint(model, coefs, set)
    return SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(
        gram_matrix, gram_constraint, zero_constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{<:AbstractAlgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(::Type{SOSPolynomialBridge{T, F, DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return [(F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT})]
end
function MOIB.concrete_bridge_type(::Type{<:SOSPolynomialBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:SOSPolynomialSet{DT, CT, <:PolyJuMP.MonomialBasis, MT, MVT}}) where {T, DT<:AbstractAlgebraicSet, CT, MT, MVT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    return SOSPolynomialBridge{T, G, DT, PolyJuMP.MonomialBasis, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::SOSPolynomialBridge{T, F, DT, BT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return 1
end
function MOI.get(b::SOSPolynomialBridge{T, F, DT, BT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::SOSPolynomialBridge)
    # First delete the constraints in which the Gram matrix appears
    MOI.delete(model, bridge.zero_constraint)
    MOI.delete(model, bridge.gram_constraint)
    # Now we delete the Gram matrix
    gram_delete(model, bridge.gram_matrix)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge)
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(model::MOI.ModelLike, ::MomentMatrix,
                 bridge::SOSPolynomialBridge)
    dual = MOI.get(model, MOI.ConstraintDual(), bridge.gram_constraint)
    monos = bridge.gram_matrix.x
    dual_matrix = MultivariateMoments.SymMatrix(dual, length(monos))
    return matmeasure(dual_matrix, monos)
end
function MOI.get(model::MOI.ModelLike, ::GramMatrix, bridge::SOSPolynomialBridge)
    return primal_value(model, bridge.gram_matrix)
end
function MOI.get(model::MOI.ModelLike, attr::CertificateMonomials,
                 bridge::SOSPolynomialBridge)
    return bridge.gram_matrix.x
end
