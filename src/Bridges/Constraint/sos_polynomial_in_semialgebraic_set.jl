"""
    SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,BT,M,NB,GB,W} <: Bridges.Constraint.AbstractBridge

`SOSPolynomialInSemialgebraicSetBridge` implements a reformulation from
[`SumOfSquares.SOSPolynomialSet{<:BasicSemialgebraicSet}`](@ref) into a
single [`SumOfSquares.WeightedSOSCone`](@ref) carrying the whole
Putinar-style decomposition
``p = \\sigma_0 + \\sum_i g_i \\sigma_i``.

For each inequality ``g_i`` of the basic semialgebraic domain, the bridge
collects a multiplier basis ``b_i`` and produces a `WeightedSOSCone` with
gram bases ``[b_0, b_1, \\dots]`` and weights ``[1, g_1, \\dots]``. The
polynomial is reduced modulo the algebraic part ``\\mathcal{V}`` before
being passed on. Whichever bridge handles `WeightedSOSCone` downstream
(`Variable.KernelBridge`, `Variable.LowRankBridge` or `Constraint.ImageBridge`)
ultimately produces the PSD constraints; this bridge adds no variables of
its own.

## Source node

`SOSPolynomialInSemialgebraicSetBridge` supports:

  * `F` in [`SumOfSquares.SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet}`](@ref)

## Target nodes

`SOSPolynomialInSemialgebraicSetBridge` creates:

  * `F` in [`SumOfSquares.WeightedSOSCone`](@ref)
"""
struct SOSPolynomialInSemialgebraicSetBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    CT<:Certificate.AbstractIdealCertificate,
    BT<:MB.SubBasis{MB.Monomial},
    M,
    NB<:SA.ExplicitBasis,
    GB<:SA.ExplicitBasis,
    W<:SA.AlgebraElement,
} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{F,SOS.WeightedSOSCone{M,NB,GB,W}}
    basis::BT
    # `multiplier_indices[1]` is the range of `gram_bases` that corresponds
    # to σ_0 (possibly multiple if the certificate uses sparsity).
    # `multiplier_indices[1+i]` is the range of `gram_bases` of the
    # Lagrangian multiplier σ_i.
    multiplier_indices::Vector{Union{Int,UnitRange{Int}}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,BT,M,NB,GB,W}},
    model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction,
    set::SOS.SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet},
) where {T,F,DT,CT,BT,M,NB,GB,W}
    @assert MOI.output_dimension(f) == length(set.basis)
    poly = MB.algebra_element(
        SA.SparseCoefficients(
            copy(collect(set.basis.keys)),
            MOI.Utilities.scalarize(f),
            SA.comparable(MB.implicit_basis(set.basis)),
        ),
        MB.implicit_basis(set.basis),
    )
    ideal_cert = Certificate.ideal_certificate(set.certificate)
    # Reduce by the algebraic ideal part `V` of the basic semialgebraic
    # domain so the rest of the certificate operates modulo it.
    domain_V = MP.similar(set.domain.V, T)
    poly_reduced = Certificate.reduced_polynomial(ideal_cert, poly, domain_V)
    sigma_0_basis = Certificate.gram_basis(
        ideal_cert,
        Certificate.with_variables(poly_reduced, set.domain),
    )
    implicit_basis = MB.implicit_basis(set.basis)
    # Build `[1, g_1, …]` weights and `[σ_0_basis, σ_1_basis, …]` gram bases
    # so the whole Putinar decomposition becomes a single `WeightedSOSCone`.
    preprocessed =
        Certificate.preprocessed_domain(set.certificate, set.domain, poly)
    # Use the same `sparse_coefficients`-based construction for the unit
    # weight as for the polynomial g_i weights so they all share the same
    # concrete element type and fit in a homogeneous `Vector{W}`.
    some_mono = first(MB.keys_as_monomials(set.basis))
    unit_poly = MP.polynomial(MP.term(one(T), MP.constant_monomial(some_mono)))
    # Append the Lagrangian multiplier bases (σ_1, σ_2, …) first and σ_0
    # last so that, downstream, `Variable.KernelBridge` allocates the σ_i
    # gram variables before the σ_0 ones — matching the variable ordering
    # of the original two-bridge pipeline (and therefore the order assumed
    # by the cached `Mock/` tests).
    gram_bases_raw = Any[]
    weights_raw = W[]
    for index in Certificate.preorder_indices(set.certificate, preprocessed)
        push!(
            gram_bases_raw,
            Certificate.multiplier_basis(set.certificate, index, preprocessed),
        )
        g = Certificate.generator(set.certificate, index, preprocessed)
        push!(
            weights_raw,
            MB.algebra_element(
                MB.sparse_coefficients(one(T) * similar(g, T)),
                implicit_basis,
            ),
        )
    end
    push!(gram_bases_raw, sigma_0_basis)
    push!(
        weights_raw,
        MB.algebra_element(MB.sparse_coefficients(unit_poly), implicit_basis),
    )
    # `gram_basis` returns either a single basis or a `Vector{basis}` when
    # sparsity is in play. Reuse `_flatten` from the SOSPolynomial bridge to
    # turn the latter into a uniform `Vector{GB}` (replicating the matching
    # weight for each sparsity block). At the same time, record where each
    # multiplier σ_i lives in the flattened `gram_bases` so that we can
    # rebuild block-diagonal gram and moment matrices per multiplier.
    # `multiplier_indices` is `[σ_0, σ_1, …]` in that order, even though the
    # underlying `gram_bases` lay them out as `[σ_1, …, σ_0]`.
    raw_ranges = Union{Int,UnitRange{Int}}[]
    cur = 0
    for raw in gram_bases_raw
        n = raw isa AbstractVector ? length(raw) : 1
        if n == 1
            cur += 1
            push!(raw_ranges, cur)
        else
            push!(raw_ranges, (cur+1):(cur+n))
            cur += n
        end
    end
    # σ_0 is the last entry of `gram_bases_raw`; bring it back to position 1
    # in `multiplier_indices` without splatting the `UnitRange` element.
    multiplier_indices = Union{Int,UnitRange{Int}}[raw_ranges[end]]
    append!(multiplier_indices, @view raw_ranges[1:(end-1)])
    gram_bases, weights, _ = _flatten(identity.(gram_bases_raw), weights_raw)
    new_basis = Certificate.zero_basis(
        ideal_cert,
        MB.explicit_basis(poly_reduced),
        domain_V,
        gram_bases,
        weights,
    )
    new_coeffs = SA.coeffs(poly_reduced, new_basis)
    constraint = MOI.add_constraint(
        model,
        MOI.Utilities.vectorize(SA.values(new_coeffs)),
        SOS.WeightedSOSCone{M}(new_basis, gram_bases, weights),
    )
    return SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,BT,M,NB,GB,W}(
        constraint,
        set.basis,
        multiplier_indices,
    )
end

function MOI.supports_constraint(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet}},
) where {T}
    return true
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:SOSPolynomialInSemialgebraicSetBridge},
)
    return Tuple{Type}[]
end
function MOI.Bridges.added_constraint_types(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,BT,M,NB,GB,W}},
) where {T,F,DT,CT,BT,M,NB,GB,W}
    return Tuple{Type,Type}[(F, SOS.WeightedSOSCone{M,NB,GB,W})]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{
        <:SOS.SOSPolynomialSet{
            SemialgebraicSets.BasicSemialgebraicSet{S,PS,AT},
            BT,
            CT,
        },
    },
) where {T,S,PS,AT,CT,BT<:MB.MonomialIndexedBasis{MB.Monomial}}
    G = MOI.Utilities.promote_operation(-, T, F, MOI.VectorOfVariables)
    IC = Certificate.ideal_certificate(CT)
    M = SOS.matrix_cone_type(IC)
    # The multiplier gram basis of the inner ideal certificate is used both
    # for σ_0 and (via `multiplier_basis_type`) for each Lagrangian σ_i.
    GB =
        SOSPolynomialBridgeType_helper_gram_basis_type(IC, MP.monomial_type(BT))
    # All weights (the σ_0 unit weight and each g_i) are built with
    # `MB.sparse_coefficients` so they share a single concrete type with
    # `Vector{Int}` keys and `Vector{T}` values.
    A = MA.promote_operation(
        MB.algebra,
        MA.promote_operation(MB.implicit_basis, BT),
    )
    C = SA.SparseCoefficients{
        Vector{Int},
        T,
        Vector{Vector{Int}},
        Vector{T},
        MP.ordering(BT),
    }
    W = SA.AlgebraElement{T,A,C}
    NB = MA.promote_operation(
        Certificate.zero_basis,
        IC,
        BT,
        SemialgebraicSets.similar_type(AT, T),
        Vector{GB},
        Vector{W},
    )
    return SOSPolynomialInSemialgebraicSetBridge{T,G,AT,IC,BT,M,NB,GB,W}
end

# A small helper to keep the `concrete_bridge_type` readable.
function SOSPolynomialBridgeType_helper_gram_basis_type(
    ::Type{IC},
    ::Type,
) where {IC}
    return _eltype(MA.promote_operation(SOS.Certificate.gram_basis, IC))
end

# Attributes, Bridge acting as a model
function MOI.get(
    ::SOSPolynomialInSemialgebraicSetBridge,
    ::MOI.NumberOfVariables,
)
    return 0
end
function MOI.get(
    ::SOSPolynomialInSemialgebraicSetBridge,
    ::MOI.ListOfVariableIndices,
)
    return MOI.VariableIndex[]
end
function MOI.get(
    ::SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,BT,M,NB,GB,W},
    ::MOI.NumberOfConstraints{F,SOS.WeightedSOSCone{M,NB,GB,W}},
) where {T,F,DT,CT,BT,M,NB,GB,W}
    return 1
end
function MOI.get(
    b::SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,BT,M,NB,GB,W},
    ::MOI.ListOfConstraintIndices{F,SOS.WeightedSOSCone{M,NB,GB,W}},
) where {T,F,DT,CT,BT,M,NB,GB,W}
    return [b.constraint]
end

# Indices
function MOI.delete(
    model::MOI.ModelLike,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    MOI.delete(model, bridge.constraint)
    return
end

# Attributes, Bridge acting as a constraint

# The monomials might be different from the ones of the original polynomial
# because of the ∑ λ_i s_i(x) so we don't define ConstraintPrimal and
# ConstraintDual, as the caller won't know how to reshape it
function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintPrimal,
    ::SOSPolynomialInSemialgebraicSetBridge,
)
    return throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    dual = MOI.get(model, attr, bridge.constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.constraint)
    μ = MultivariateMoments.moment_vector(dual, set.basis)
    monos = MB.keys_as_monomials(bridge.basis)
    return [dot(mono, μ) for mono in monos]
end
function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    # We can't forward `MomentsAttribute` to `bridge.constraint` because
    # nothing in the `WeightedSOSCone` chain knows about
    # `MomentsAttribute`. Build it directly from the dual and the new basis
    # stored on the `WeightedSOSCone` set.
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
    attr::SOS.CertificateBasis,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    return MOI.get(model, attr, bridge.constraint)
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{SOS.GramMatrixAttribute,SOS.MomentMatrixAttribute},
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    # `multiplier_index = 0` means "σ_0"; remap it to whichever indices the
    # σ_0 block(s) occupy in the underlying `WeightedSOSCone`'s
    # `gram_bases`.
    SOS.check_multiplier_index_bounds(attr, 0:0)
    return _get(model, attr, bridge.constraint, bridge.multiplier_indices[1])
end

function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.LagrangianMultipliers,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    # `multiplier_indices[1]` is σ_0; the rest are σ_1, σ_2, … in order.
    # Each multiplier is returned as a `BlockDiagonalGramMatrix`, even when
    # its `gram_bases` range has a single entry, so that the result type is
    # uniform across multipliers.
    return map(2:length(bridge.multiplier_indices)) do i
        idx = bridge.multiplier_indices[i]
        range = idx isa Int ? (idx:idx) : idx
        gram_attr = SOS.GramMatrixAttribute(; result_index = attr.result_index)
        return _get(model, gram_attr, bridge.constraint, range)
    end
end
