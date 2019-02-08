function moi_matpoly(model::MOI.ModelLike, monos)
    return  MatPolynomial{MOI.SingleVariable}((i, j) -> MOI.SingleVariable(MOI.add_variable(model)),
                                              monos)
end
function matpoly_in_cone(model::MOI.ModelLike, monos, set::SOSLikeCones)
    p = moi_matpoly(model, monos)
    matrix_add_constraint(model, p, matrix_cone(set))
    return p
end
function _matposynomial(model::MOI.ModelLike, monos)
    p = moi_matpoly(model, monos)
    for q in p.Q.Q
        MOI.add_constraint(model, q, MOI.GreaterThan(0.0))
    end
    return p
end
function matpoly_in_cone(model, x, set::CoSOSLikeCones)
    _matplus(matpoly_in_cone(model, x, _nococone(set)), _matposynomial(m, x))
end

struct SOSPolynomialBridge{T, F <: MOI.AbstractVectorFunction,
                           DT <: AbstractSemialgebraicSet,
                           BT <: PolyJuMP.AbstractPolynomialBasis,
                           MT <: AbstractMonomial,
                           MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    gram_matrix::MatPolynomial{MOI.SingleVariable, MT, MVT}
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}
end

function SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(model::MOI.ModelLike,
                                                    f::MOI.AbstractVectorFunction,
                                                    s::SOSPolynomialSet{<:AbstractAlgebraicSet}) where {T, F, DT, BT, MT, MVT}
    @assert MOI.output_dimension(f) == length(s.monomials)
    p = polynomial(collect(MOIU.eachscalar(f)), s.monomials)
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = convert(typeof(p), rem(p, ideal(s.domain)))
    X = monomials_half_newton_polytope(monomials(r), s.newton_polytope)
    gram_matrix = matpoly_in_cone(model, X, s.cone)
    q = r - gram_matrix
    set = PolyJuMP.ZeroPolynomialSet(s.domain, s.basis, monomials(q))
    zero_constraint = MOI.add_constraint(model, MOIU.vectorize(coefficients(q)),
                                         set)
    return SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(gram_matrix, zero_constraint)
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
    # Now we delete the Gram matrix
    for vi in bridge.gram_matrix.Q.Q
        MOI.delete(model, vi.variables)
    end
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge)
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(model::MOI.ModelLike, ::MomentMatrix,
                 bridge::SOSPolynomialBridge)
    return primal_value(model, bridge.gram_matrix)
    μ = MOI.get(model, MOI.ConstraintDual(), bridge)
end
function MOI.get(model::MOI.ModelLike, ::GramMatrix, bridge::SOSPolynomialBridge)
    return primal_value(model, bridge.gram_matrix)
end
function MOI.get(model::MOI.ModelLike, attr::CertificateMonomials,
                 bridge::SOSPolynomialBridge)
    return bridge.gram_matrix.x
end
