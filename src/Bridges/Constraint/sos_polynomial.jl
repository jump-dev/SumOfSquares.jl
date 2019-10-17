struct SOSPolynomialBridge{
    T, F <: MOI.AbstractVectorFunction,
    DT <: SemialgebraicSets.AbstractSemialgebraicSet,
    UMCT <: MOI.ConstraintIndex{MOI.VectorOfVariables},
    UMST, MCT,
    BT <: PolyJuMP.AbstractPolynomialBasis,
    CT <: SOS.Certificate.AbstractIdealCertificate,
    MT <: MP.AbstractMonomial,
    MVT <: AbstractVector{MT}} <: MOIB.Constraint.AbstractBridge

    Q::Vector{MOI.VariableIndex} # TODO the type will be different for sparse SOS
    cQ::UMCT # TODO the type will be different for sparse SOS
    certificate_monomials::MVT
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}
    domain::DT
    monomials::MVT
    certificate::CT
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT}},
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}) where {T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT}

    @assert MOI.output_dimension(f) == length(s.monomials)
    p = MP.polynomial(collect(MOIU.eachscalar(f)), s.monomials)
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
    # need to call `changecoefficienttype`. This is critical since `T` is
    # `Float64` when used with JuMP and the coefficient type is often `Int` if
    # `set.domain.V` is `FullSpace` or `FixedPolynomialsSet`.
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = SOS.Certificate.get(s.certificate, SOS.Certificate.ReducedPolynomial(), p, MP.changecoefficienttype(s.domain, T))
    X = SOS.Certificate.get(s.certificate, SOS.Certificate.GramBasis(), r)
    g, Q, cQ = SOS.add_gram_matrix(model, MCT, X)
    q = r - g
    set = PolyJuMP.ZeroPolynomialSet(s.domain, SOS.Certificate.zero_basis(s.certificate), MP.monomials(q))
    coefs = MOIU.vectorize(MP.coefficients(q))
    zero_constraint = MOI.add_constraint(model, coefs, set)
    return SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT}(
        Q, cQ, X, zero_constraint, s.domain, s.monomials, s.certificate)
end

function MOI.supports_constraint(::Type{SOSPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}}) where T
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT}}) where {T, F, DT, UMCT, UMST, MCT}
    return constrained_variable_types(MCT)
end
function MOIB.added_constraint_types(::Type{SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT}}) where {T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT}
    return [(F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT})]
end
function MOIB.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{DT, MT, MVT, CT}}) where {T, DT<:SemialgebraicSets.AbstractAlgebraicSet, MT, MVT, CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    MCT = SOS.matrix_cone_type(CT)
    UMCT = union_constraint_types(MCT)
    UMST = union_set_types(MCT)
    BT = SOS.Certificate.zero_basis_type(CT)
    return SOSPolynomialBridge{T, G, DT, UMCT, UMST, MCT, BT, CT, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::SOSPolynomialBridge, ::MOI.NumberOfVariables)
    return length(bridge.Q)
end
function MOI.get(bridge::SOSPolynomialBridge, ::MOI.ListOfVariableIndices)
    return bridge.Q
end
function MOI.get(bridge::SOSPolynomialBridge{T, F, DT, UMCT, UMST},
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, S}) where {T, F, DT, UMCT, UMST, S<:UMST}
    return bridge.cQ isa MOI.ConstraintIndex{MOI.VectorOfVariables, S} ? 1 : 0
end
function MOI.get(bridge::SOSPolynomialBridge{T, F, DT, UMCT, UMST},
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}) where {T, F, DT, UMCT, UMST, S<:UMST}
    if bridge.cQ isa MOI.ConstraintIndex{MOI.VectorOfVariables, S}
        return [bridge.cQ]
    else
        return MOI.ConstraintIndex{MOI.VectorOfVariables, S}[]
    end
end
function MOI.get(::SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {
        T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT
    }
    return 1
end
function MOI.get(b::SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {
        T, F, DT, UMCT, UMST, MCT, BT, CT, MT, MVT
    }
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::SOSPolynomialBridge)
    # First delete the constraints in which the Gram matrix appears
    MOI.delete(model, bridge.zero_constraint)
    # Now we delete the Gram matrix
    MOI.delete(model, bridge.Q)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 ::MOI.ConstraintSet,
                 bridge::SOSPolynomialBridge)
    return SOS.SOSPolynomialSet(bridge.domain, bridge.monomials, bridge.certificate)
end
function MOI.get(::MOI.ModelLike,
                 ::MOI.ConstraintPrimal,
                 ::SOSPolynomialBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge)
    dual = MOI.get(model, attr, bridge.zero_constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.zero_constraint)
    μ = MultivariateMoments.measure(dual, set.monomials)
    I = SemialgebraicSets.ideal(bridge.domain)
    return [dot(rem(mono, I), μ) for mono in bridge.monomials]
end
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge{T, <:MOI.AbstractVectorFunction, SemialgebraicSets.FullSpace}) where {T}
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(model::MOI.ModelLike, attr::PolyJuMP.MomentsAttribute,
                 bridge::SOSPolynomialBridge)
    return MOI.get(model, attr, bridge.zero_constraint)
end

function MOI.get(::MOI.ModelLike, ::SOS.CertificateMonomials,
                 bridge::SOSPolynomialBridge)
    return bridge.certificate_monomials
end
function MOI.get(model::MOI.ModelLike,
                 attr::SOS.GramMatrixAttribute,
                 bridge::SOSPolynomialBridge)
    return SOS.build_gram_matrix(MOI.get(model, MOI.VariablePrimal(attr.N), bridge.Q),
                                 bridge.certificate_monomials)
end
function MOI.get(model::MOI.ModelLike,
                 attr::SOS.MomentMatrixAttribute,
                 bridge::SOSPolynomialBridge)
    return SOS.build_moment_matrix(MOI.get(model, MOI.ConstraintDual(attr.N),
                                           bridge.cQ),
                                   bridge.certificate_monomials)
end
