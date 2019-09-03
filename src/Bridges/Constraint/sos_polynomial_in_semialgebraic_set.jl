function lagrangian_multiplier(model::MOI.ModelLike, p, certificate, index, domain, T::Type)
    monos = Certificate.get(certificate, Certificate.MultiplierBasis(), index, domain, p)
    MCT = SOS.matrix_cone_type(typeof(certificate))
    Q, con_Q = MOI.add_constrained_variables(model, SOS.matrix_cone(MCT, length(monos)))
    return SOS.build_gram_matrix(Q, monos), Q, con_Q, monos
end

struct SOSPolynomialInSemialgebraicSetBridge{
    T, F <: MOI.AbstractVectorFunction, DT <: SemialgebraicSets.AbstractSemialgebraicSet,
    CT <: Certificate.AbstractIdealCertificate,
    UMCT <: MOI.ConstraintIndex{MOI.VectorOfVariables},
    UMST,
    MT <: MP.AbstractMonomial, MVT <: AbstractVector{MT}} <: MOIB.Constraint.AbstractBridge
    lagrangian_monomials::Vector{MVT}
    lagrangian_variables::Vector{Vector{MOI.VariableIndex}}
    lagrangian_constraints::Vector{UMCT}
    constraint::MOI.ConstraintIndex{F, SOS.SOSPolynomialSet{DT, MT, MVT, CT}}
    monomials::MVT
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, UMCT, UMST, MT, MVT}},
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    set::SOS.SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet}) where {T, F, DT, CT, UMCT, UMST, MT, MVT}

    @assert MOI.output_dimension(f) == length(set.monomials)
    p = MP.polynomial(collect(MOIU.eachscalar(f)), set.monomials)
    n = length(set.domain.p)
    λ_monos     = MVT[]
    λ_variables = Vector{MOI.VariableIndex}[]
    λ_constraints = UMCT[]
    for index in Certificate.get(set.certificate, Certificate.PreorderIndices(), set.domain)
        λ, λ_variable, λ_constraint, λ_mono = lagrangian_multiplier(
            model, p, set.certificate, index, set.domain, T)
        push!(λ_variables, λ_variable)
        push!(λ_constraints, λ_constraint)
        push!(λ_monos, λ_mono)
        # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
        # need to call `changecoefficienttype`. This is critical since `T` is
        # `Float64` when used with JuMP and the coefficient type is often `Int` if
        # `set.domain.V` is `FullSpace` or `FixedPolynomialsSet`.
        g = Certificate.get(set.certificate, Certificate.Generator(), index, set.domain)
        p -= λ * MP.changecoefficienttype(g, T)
    end
    new_set = SOS.SOSPolynomialSet(
        set.domain.V, MP.monomials(p),
        Certificate.get(set.certificate, Certificate.IdealCertificate()))
    constraint = MOI.add_constraint(model, MOIU.vectorize(MP.coefficients(p)),
                                    new_set)

    return SOSPolynomialInSemialgebraicSetBridge{
        T, F, DT, CT, UMCT, UMST, MT, MVT}(λ_monos, λ_variables, λ_constraints, constraint, set.monomials)
end

function MOI.supports_constraint(::Type{SOSPolynomialInSemialgebraicSetBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOS.SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet}}) where T
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT}}) where {T, F, DT, CT}
    return constrained_variable_types(SOS.matrix_cone_type(CT))
end
function MOIB.added_constraint_types(
    ::Type{SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, UMCT, UMST, MT, MVT}}) where {T, F, DT, CT, UMCT, UMST, MT, MVT}
    return [(F, SOS.SOSPolynomialSet{DT, MT, MVT, CT})]
end
function MOIB.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{SemialgebraicSets.BasicSemialgebraicSet{S, PS, AT}, MT, MVT, CT}}) where {T, S, PS, AT, CT, MT, MVT}

    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    MCT = SOS.matrix_cone_type(CT)
    UMCT = union_constraint_types(MCT)
    UMST = union_set_types(MCT)
    IC = Certificate.get(CT, Certificate.IdealCertificate())
    return SOSPolynomialInSemialgebraicSetBridge{T, G, AT, IC, UMCT, UMST, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::SOSPolynomialInSemialgebraicSetBridge,
                 ::MOI.NumberOfVariables)
    return reduce((n, vis) -> n + length(vis), bridge.lagrangian_variables, init=0)
end
function MOI.get(bridge::SOSPolynomialInSemialgebraicSetBridge, ::MOI.ListOfVariableIndices)
    return Iterators.flatten(bridge.lagrangian_variables)
end
function MOI.get(bridge::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, UMCT, UMST},
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, S}) where {T, F, DT, CT, UMCT, UMST, S<:UMST}
    return reduce((n, ci) -> n + (ci isa MOI.ConstraintIndex{MOI.VectorOfVariables, S}),
                  bridge.lagrangian_constraints, init=0)
end
function MOI.get(b::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, UMCT, UMST},
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}) where {T, F, DT, CT, UMCT, UMST, S<:UMST}
    return filter(ci -> ci isa MOI.ConstraintIndex{MOI.VectorOfVariables, S}, bridge.lagrangian_constraints)
end
function MOI.get(::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, UMCT, UMST, MT, MVT},
                 ::MOI.NumberOfConstraints{F, SOS.SOSPolynomialSet{DT, MT, MVT, CT}}) where {T, F, DT, CT, UMCT, UMST, MT, MVT}
    return 1
end
function MOI.get(b::SOSPolynomialInSemialgebraicSetBridge{T, F, DT, CT, UMCT, UMST, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, SOS.SOSPolynomialSet{DT, MT, MVT, CT}}) where {T, F, DT, CT, UMCT, UMST, MT, MVT}
    return [b.constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike,
                    bridge::SOSPolynomialInSemialgebraicSetBridge)
    MOI.delete(model, bridge.constraint)
    for variables in bridge.lagrangian_variables
        MOI.delete(model, variables)
    end
end

# Attributes, Bridge acting as a constraint

# The monomials might be different from the ones of the original polynomial
# because of the ∑ λ_i s_i(x) so we don't define ConstraintPrimal and
# ConstraintDual, as the caller won't know how to reshape it
function MOI.get(::MOI.ModelLike,
                 ::MOI.ConstraintPrimal,
                 ::SOSPolynomialInSemialgebraicSetBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    dual = MOI.get(model, attr, bridge.constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.constraint)
    μ = MultivariateMoments.measure(dual, set.monomials)
    return [dot(mono, μ) for mono in bridge.monomials]
end
function MOI.get(model::MOI.ModelLike, attr::PolyJuMP.MomentsAttribute,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, bridge.constraint)
end

function MOI.get(model::MOI.ModelLike,
                 attr::Union{SOS.CertificateMonomials, SOS.GramMatrixAttribute,
                             SOS.MomentMatrixAttribute},
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    return MOI.get(model, attr, bridge.constraint)
end
function MOI.get(model::MOI.ModelLike, attr::SOS.LagrangianMultipliers,
                 bridge::SOSPolynomialInSemialgebraicSetBridge)
    @assert eachindex(bridge.lagrangian_variables) == eachindex(bridge.lagrangian_monomials)
    map(i -> SOS.build_gram_matrix(MOI.get(model, MOI.VariablePrimal(attr.N),
                                           bridge.lagrangian_variables[i]),
                                   bridge.lagrangian_monomials[i]),
        eachindex(bridge.lagrangian_variables))
end
