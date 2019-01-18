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
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}
end

function SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(model::MOI.ModelLike,
                                                    f::MOI.AbstractVectorFunction,
                                                    s::SOSPolynomialSet{<:AbstractAlgebraicSet}) where {T, F, DT, BT, MT, MVT}
    @assert MOI.output_dimension(f) == length(s.monomials)
    p = polynomial(collect(MOIU.eachscalar(f)), s.monomials)
    r = rem(p, ideal(s.domain))
    X = monomials_half_newton_polytope(monomials(r), s.newton_polytope)
    slack = matpoly_in_cone(model, X, s.cone)
    q = r - slack
    set = PolyJuMP.ZeroPolynomialSet(s.domain, s.basis, monomials(q))
    zero_constraint = MOI.add_constraint(model, MOIU.vectorize(coefficients(q)),
                                         set)
    return SOSPolynomialBridge{T, F, DT, BT, MT, MVT}(zero_constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{FullSpace}}) where T
    return true
end
function MOIB.added_constraint_types(::Type{SOSPolynomialBridge{T, F, DT, BT, MT, MVT}}) where {T, F, DT, BT, MT, MVT}
    return [(F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT})]
end
function MOIB.concrete_bridge_type(::Type{<:SOSPolynomialBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:SOSPolynomialSet{FullSpace, CT, <:PolyJuMP.MonomialBasis, MT, MVT}}) where {T, CT, MT, MVT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    return SOSPolynomialBridge{T, G, FullSpace, PolyJuMP.MonomialBasis, MT, MVT}
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
