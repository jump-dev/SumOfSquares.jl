struct SOSPolynomialBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    M,
    G<:MB.AbstractPolynomialBasis,
    CT<:SOS.Certificate.AbstractIdealCertificate,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
    W<:MP.AbstractTerm{T},
} <: MOI.Bridges.Constraint.SetMapBridge{T,SOS.WeightedSOSCone{M,MB.SubBasis{MB.Monomial,MT,MVT},G,W},SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT},F,F}
    constraint::MOI.ConstraintIndex{F,SOS.WeightedSOSCone{M,MB.SubBasis{MB.Monomial,MT,MVT},G,W}}
    set::SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialBridge{T,F,DT,M,G,CT,MT,MVT,W}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet},
) where {T,F,DT,M,G,CT,MT,MVT,W}
    @assert MOI.output_dimension(func) == length(set.monomials)
    # MOI does not modify the coefficients of the functions so we can modify `p`.
    # without altering `f`.
    # The monomials may be copied by MA however so we need to copy it.
    p = MP.polynomial(MOI.Utilities.scalarize(func), copy(set.monomials))
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
    # need to call `similar`. This is critical since `T` is
    # `Float64` when used with JuMP and the coefficient type is often `Int` if
    # `set.domain.V` is `FullSpace` or `FixedPolynomialSet`.
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = SOS.Certificate.reduced_polynomial(
        set.certificate,
        p,
        MP.similar(set.domain, T),
    )
    gram_basis = SOS.Certificate.gram_basis(
        set.certificate,
        SOS.Certificate.with_variables(r, set.domain),
    )
    constraint = MOI.add_constraint(
        model,
        func,
        SOS.WeightedSOSCone{M}(
            MB.SubBasis{MB.Monomial}(set.monomials),
            [gram_basis],
            [MP.term(one(T), MP.constant_monomial(p))],
        ),
    )
    return SOSPolynomialBridge{T,F,DT,M,G,CT,MT,MVT,W}(
        constraint,
        set,
    )
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
    G = SOS.Certificate.gram_basis_type(CT, MT)
    W = MP.term_type(MT, T)
    return SOSPolynomialBridge{T,F,DT,M,G,CT,MT,MVT,W}
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
    return MultivariateMoments.Measure(
        MOI.get(model, MOI.ConstraintDual(attr.result_index), bridge.constraint),
        bridge.set.monomials,
    )
end

function MOI.get(
    ::MOI.ModelLike,
    ::SOS.CertificateBasis,
    bridge::SOSPolynomialBridge,
)
    return bridge.gram_basis
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{SOS.GramMatrixAttribute,SOS.MomentMatrixAttribute,SOS.SOSDecompositionAttribute},
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
