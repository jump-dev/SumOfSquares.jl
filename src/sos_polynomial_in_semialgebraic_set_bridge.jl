function lagrangian_multiplier(model::MOI.ModelLike, p, certificate, index, domain, T::Type)
    monos = Certificate.get(certificate, Certificate.MultiplierBasis(), index, domain, p)
    Q, variable_bridge = add_matrix_variable_bridge(
        model, matrix_cone_type(typeof(certificate)), length(monos), T)
    return build_gram_matrix(Q, monos), variable_bridge, monos
end

struct SOSPolynomialInSemialgebraicSetBridge{
    T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
    CT <: Certificate.AbstractIdealCertificate, VBS <: AbstractVariableBridge,
    MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    lagrangian_monomials::Vector{MVT}
    lagrangian_bridges::Vector{VBS}
    constraint::MOI.ConstraintIndex{F, SOSPolynomialSet{DT, MT, MVT, CT}}
    monomials::MVT
end

function SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, VBS, MT, MVT}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    set::SOSPolynomialSet{<:BasicSemialgebraicSet}) where {T, F, DT, CT, VBS, MT, MVT}
    @assert MOI.output_dimension(f) == length(set.monomials)
    p = MP.polynomial(collect(MOIU.eachscalar(f)), set.monomials)
    n = length(set.domain.p)
    λ_monos   = MVT[]
    λ_bridges = VBS[]
    for index in Certificate.get(set.certificate, Certificate.PreorderIndices(), set.domain)
        λ, λ_bridge, λ_mono = lagrangian_multiplier(
            model, p, set.certificate, index, set.domain, T)
        push!(λ_bridges, λ_bridge)
        push!(λ_monos, λ_mono)
        # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
        # need to call `changecoefficienttype`. This is critical since `T` is
        # `Float64` when used with JuMP and the coefficient type is often `Int` if
        # `set.domain.V` is `FullSpace` or `FixedPolynomialsSet`.
        g = Certificate.get(set.certificate, Certificate.Generator(), index, set.domain)
        p -= λ * MP.changecoefficienttype(g, T)
    end
    new_set = SOSPolynomialSet(
        set.domain.V, MP.monomials(p),
        Certificate.get(set.certificate, Certificate.IdealCertificate()))
    constraint = MOI.add_constraint(model, MOIU.vectorize(MP.coefficients(p)),
                                    new_set)

    return SOSPolynomialInSemialgebraicSetBridge{
        T, F, DT, CT, VBS, MT, MVT}(λ_monos, λ_bridges, constraint, set.monomials)
end

function MOI.supports_constraint(::Type{SOSPolynomialInSemialgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{<:BasicSemialgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, VBS, MT, MVT}}) where {T, F, DT, CT, VBS, MT, MVT}
    added = [(F, SOSPolynomialSet{DT, MT, MVT, CT})]
    return append_added_constraint_types(added, matrix_cone_type(CT), T)
end
function MOIB.concrete_bridge_type(::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:SOSPolynomialSet{BasicSemialgebraicSet{S, PS, AT}, MT, MVT, CT}}) where {T, S, PS, AT, CT, MT, MVT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    VBS = union_vector_bridge_types(matrix_cone_type(CT), T)
    IC = Certificate.get(CT, Certificate.IdealCertificate())
    return SOSPolynomialInSemialgebraicSetBridge{T, G, AT, IC, VBS, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, VBS, MT, MVT},
                 ::MOI.NumberOfConstraints{F, SOSPolynomialSet{DT, MT, MVT, CT}}) where {T, F, DT, CT, VBS, MT, MVT}
    return 1
end
function MOI.get(b::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, VBS, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, SOSPolynomialSet{DT, MT, MVT, CT}}) where {T, F, DT, CT, VBS, MT, MVT}
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
function MOI.get(::MOI.ModelLike,
                 ::MOI.ConstraintPrimal,
                 ::SOSPolynomialInSemialgebraicSetBridge)
    throw(ValueNotSupported())
end

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    dual = MOI.get(model, attr, bridge.constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.constraint)
    μ = measure(dual, set.monomials)
    return [dot(mono, μ) for mono in bridge.monomials]
end
function MOI.get(model::MOI.ModelLike, attr::PolyJuMP.MomentsAttribute,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, bridge.constraint)
end

function MOI.get(model::MOI.ModelLike,
                 attr::Union{CertificateMonomials, GramMatrixAttribute,
                             MomentMatrixAttribute},
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, bridge.constraint)
end
function MOI.get(model::MOI.ModelLike, ::LagrangianMultipliers,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    @assert eachindex(bridge.lagrangian_bridges) == eachindex(bridge.lagrangian_monomials)
    map(i -> build_gram_matrix(MOI.get(model, MOI.ConstraintPrimal(),
                                       bridge.lagrangian_bridges[i]),
                               bridge.lagrangian_monomials[i]),
        eachindex(bridge.lagrangian_bridges))
end
