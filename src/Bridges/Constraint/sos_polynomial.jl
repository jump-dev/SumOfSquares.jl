struct SOSPolynomialBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    M, # matrix cone type
    BT,
    B,
    G<:SA.ExplicitBasis,
    CT<:SOS.Certificate.AbstractIdealCertificate,
    W<:SA.AlgebraElement,
} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    SOS.WeightedSOSCone{M,B,G,W},
    SOS.SOSPolynomialSet{DT,BT,CT},
    F,
    F,
}
    constraint::MOI.ConstraintIndex{F,SOS.WeightedSOSCone{M,B,G,W}}
    set::SOS.SOSPolynomialSet{DT,BT,CT}
    new_basis::B
    flat_indices::Union{Int,Base.UnitRange{Int}}
end

function _flatten(gram_bases::Vector{<:SA.AbstractBasis}, weights)
    return gram_bases, weights, 1
end

function _flatten(
    gram_bases::Vector{Vector{B}},
    weights,
) where {B<:SA.AbstractBasis}
    flat_gram_bases = eltype(eltype(gram_bases))[]
    flat_weights = eltype(weights)[]
    for (g, w) in zip(gram_bases, weights)
        for flat in g
            push!(flat_gram_bases, flat)
            push!(flat_weights, w)
        end
    end
    return flat_gram_bases, flat_weights, 1:length(flat_weights)
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialBridge{T,F,DT,M,BT,B,G,CT,W}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet},
) where {T,F,DT,M,BT,B,G,CT,W}
    @show func
    @show set.basis
    @assert MOI.output_dimension(func) == length(set.basis)
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
    # need to call `similar`. This is critical since `T` is
    # `Float64` when used with JuMP and the coefficient type is often `Int` if
    # `set.domain.V` is `FullSpace` or `FixedPolynomialSet`.
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    domain = MP.similar(set.domain, T)
    poly = SOS.Certificate.reduced_polynomial(
        set.certificate,
        # MOI does not modify the coefficients of the functions so we can modify `p`.
        # without altering `f`.
        # The basis may be copied by MA however so we need to copy it.
        MB.algebra_element(MOI.Utilities.scalarize(func), copy(set.basis)),
        domain,
    )
    gram_basis = SOS.Certificate.gram_basis(
        set.certificate,
        SOS.Certificate.with_variables(poly, set.domain),
    )
    gram_bases = [gram_basis]
    weights = [MB.constant_algebra_element(typeof(SA.basis(poly)), T)]
    @show typeof(gram_bases)
    flat_gram_bases, flat_weights, flat_indices = _flatten(gram_bases, weights)
    @show flat_indices
    new_basis = SOS.Certificate.reduced_basis(
        set.certificate,
        SA.basis(poly),
        domain,
        flat_gram_bases,
        flat_weights,
    )
    new_coeffs = SA.coeffs(poly, new_basis)
    @show new_coeffs
    constraint = MOI.add_constraint(
        model,
        MOI.Utilities.vectorize(new_coeffs),
        SOS.WeightedSOSCone{M}(new_basis, flat_gram_bases, flat_weights),
    )
    return SOSPolynomialBridge{T,F,DT,M,BT,B,G,CT,W}(
        constraint,
        set,
        new_basis,
        flat_indices,
    )
end

function MOI.supports_constraint(
    ::Type{SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}},
) where {T}
    return MOI.Utilities.is_coefficient_type(F, T)
end

_eltype(::Type{Vector{T}}) where T = T
_eltype(::Type{T}) where T = T

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{DT,BT,CT}},
) where {T,DT<:SemialgebraicSets.AbstractAlgebraicSet,BT,CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    M = SOS.matrix_cone_type(CT)
    W = SA.AlgebraElement{MB.algebra_type(BT),T,Vector{T}}
    B = MA.promote_operation(
        SOS.Certificate.reduced_basis,
        CT,
        BT,
        SemialgebraicSets.similar_type(DT, T),
        Vector{BT},
        Vector{W},
    )
    G = SOS.Certificate.gram_basis_type(CT)
    return SOSPolynomialBridge{T,F,DT,M,BT,B,_eltype(G),CT,W}
end

function MOI.Bridges.inverse_map_function(::SOSPolynomialBridge, f)
    throw(MOI.Bridges.MapNotInvertible())
    # Does not work with QuotientBasis
    #return SA.coeffs(MP.polynomial(f, bridge.new_basis), bridge.set.basis)
end

function MOI.Bridges.adjoint_map_function(bridge::SOSPolynomialBridge, f)
    # FIXME `coeffs` should be an `AbstractMatrix`
    return MB.adjoint_coeffs(f, bridge.new_basis, bridge.set.basis)
end

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

function _get(
    model::MOI.ModelLike,
    attr,
    constraint::MOI.ConstraintIndex,
    index::Int,
)
    return MOI.get(
        model,
        typeof(attr)(
            multiplier_index = index,
            result_index = attr.result_index,
        ),
        constraint,
    )
end

function _get(
    model::MOI.ModelLike,
    attr,
    constraint::MOI.ConstraintIndex,
    indices::UnitRange,
)
    return MultivariateMoments.block_diagonal([
        _get(model, attr, constraint, index) for index in indices
    ])
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
    return _get(model, attr, bridge.constraint, bridge.flat_indices)
end
