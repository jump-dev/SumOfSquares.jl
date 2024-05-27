struct SOSPolynomialBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    M, # matrix cone type
    B,
    G<:SA.ExplicitBasis,
    CT<:SOS.Certificate.AbstractIdealCertificate,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
    W<:MP.AbstractTerm{T},
} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    SOS.WeightedSOSCone{M,B,G,W},
    SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT},
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,SOS.WeightedSOSCone{M,B,G,W}}
    set::SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialBridge{T,F,DT,M,B,G,CT,MT,MVT,W}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet},
) where {T,F,DT,M,B,G,CT,MT,MVT,W}
    @assert MOI.output_dimension(func) == length(set.basis)
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
    # need to call `similar`. This is critical since `T` is
    # `Float64` when used with JuMP and the coefficient type is often `Int` if
    # `set.domain.V` is `FullSpace` or `FixedPolynomialSet`.
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    domain = MP.similar(set.domain, T)
    coeffs, basis = SOS.Certificate.reduced_polynomial(set.certificate, func, set.basis, domain)
    gram_basis = SOS.Certificate.gram_basis(
        set.certificate,
        SOS.Certificate.with_variables(basis, set.domain),
    )
    gram_bases = [gram_basis]
    weights = [MP.term(one(T), MP.constant_monomial(eltype(basis.monomials)))]
    constraint = MOI.add_constraint(
        model,
        coeffs,
        SOS.WeightedSOSCone{M}(
            SOS.Certificate.reduced_basis(set.certificate, basis, domain, gram_bases, weights),
            gram_bases,
            weights,
        ),
    )
    return SOSPolynomialBridge{T,F,DT,M,B,G,CT,MT,MVT,W}(constraint, set)
end

function MOI.supports_constraint(
    ::Type{SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}},
) where {T}
    return MOI.Utilities.is_coefficient_type(F, T)
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT}},
) where {T,DT<:SemialgebraicSets.AbstractAlgebraicSet,MT,MVT,CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    M = SOS.matrix_cone_type(CT)
    B = MA.promote_operation(
        SOS.Certificate.reduced_basis,
        CT,
        MB.SubBasis{MB.Monomial,MT,MVT},
        SemialgebraicSets.similar_type(DT, T),
        Vector{MB.SubBasis{MB.Monomial,MT,MVT}},
        Vector{MP.term_type(MT, T)},
    )
    G = SOS.Certificate.gram_basis_type(CT, MT)
    W = MP.term_type(MT, T)
    return SOSPolynomialBridge{T,F,DT,M,B,G,CT,MT,MVT,W}
end

MOI.Bridges.inverse_map_function(::Type{<:SOSPolynomialBridge}, f) = f
MOI.Bridges.adjoint_map_function(::Type{<:SOSPolynomialBridge}, f) = f

# Attributes, Bridge acting as a constraint
function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::SOSPolynomialBridge,
)
    return bridge.set
end

function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::SOSPolynomialBridge,
)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.constraint)
    return MultivariateMoments.moment_vector(
        MOI.get(
            model,
            MOI.ConstraintDual(attr.result_index),
            bridge.constraint,
        ),
        set.basis,
    )
end

function MOI.get(
    model::MOI.ModelLike,
    ::SOS.CertificateBasis,
    bridge::SOSPolynomialBridge,
)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.constraint)
    return set.gram_bases[]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{
        SOS.GramMatrixAttribute,
        SOS.MomentMatrixAttribute,
        SOS.SOSDecompositionAttribute,
    },
    bridge::SOSPolynomialBridge,
)
    SOS.check_multiplier_index_bounds(attr, 0:0)
    return MOI.get(
        model,
        typeof(attr)(
            multiplier_index = attr.multiplier_index + 1,
            result_index = attr.result_index,
        ),
        bridge.constraint,
    )
end
