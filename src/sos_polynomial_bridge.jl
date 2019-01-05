function _createslack(model, x, set::SOSLikeCones)
    JuMP.add_variable(model,
                      PolyJuMP.Variable(_varconetype(set)(x), false, false))
end
function _matposynomial(m, x)
    p = _matpolynomial(m, x, false, false)
    for q in p.Q
        JuMP.set_lower_bound(q, 0)
    end
    p
end
function _createslack(m, x, set::CoSOSLikeCones)
    _matplus(_createslack(m, x, _nococone(set)), _matposynomial(m, x))
end

struct SOSPolynomialBridge{T, F <: MOI.AbstractVectorFunction,
                           DT <: AbstractSemialgebraicSet,
                           BT <: PolyJuMP.AbstractPolynomialBasis,
                           MT <: AbstractMonomial,
                           MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}
end

function SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(model::MOI.ModelLike,
                                                    f::MOI.AbstractVectorFunction,
                                                    s::SOSPolynomialSet{<:AbstractAlgebraicSet}) where {T, F, DT, BT, MT, MVT}
    @assert MOI.output_dimension(f) == length(s.monomials)
    p = polynomial(MOIU.eachscalar(f), s.monomials)
    r = rem(p, ideal(domain))
    X = monomials_half_newton_polytope(monomials(r), s.newton_polytope)
    slack = _createslack(model, X, set)
    q = r - slack
    set = PolyJuMP.ZeroPolynomialSet(s.domain, s.basis, monomials(q))
    zero_constraint = MOI.add_constraint(model, coefficients(q), set)
    return SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(zero_constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{FullSpace}}) where T
    return true
end
function added_constraint_types(::Type{SOSPolynomialBridge{T, F, DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return [(F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT})]
end
function concrete_bridge_type(::Type{<:SOSPolynomialBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              ::Type{SOSPolynomialSet{FullSpace, <:MonomialBasis, MT, MVT}}) where {T, MT, MVT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, zeros(T))
    return SOSPolynomialBridge{T, G, MT, MVT}
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
function MOI.delete(model::MOI.ModelLike, c::SOSPolynomialBridge)
    MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 attr::Union{MOI.ConstraintPrimal, MOI.ConstraintDual},
                 c::SOSPolynomialBridge)
    return MOI.get(model, attr, c.zero_constraint)
end
