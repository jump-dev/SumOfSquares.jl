function lagrangian_multiplier(
    model::MOI.ModelLike,
    certificate,
    index,
    preprocessed,
    T::Type,
)
    basis = Certificate.multiplier_basis(certificate, index, preprocessed)
    MCT = SOS.matrix_cone_type(typeof(certificate))
    return SOS.add_gram_matrix(model, MCT, basis, T)..., basis
end

struct SOSPolynomialInSemialgebraicSetBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    CT<:Certificate.AbstractIdealCertificate,
    B<:Union{Vector{<:MB.SA.ExplicitBasis},MB.SA.ExplicitBasis},
    UMCT<:Union{
        Vector{<:MOI.ConstraintIndex{MOI.VectorOfVariables}},
        MOI.ConstraintIndex{MOI.VectorOfVariables},
    },
    UMST,
    MCT,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
} <: MOI.Bridges.Constraint.AbstractBridge
    lagrangian_bases::Vector{B}
    lagrangian_variables::Vector{
        Union{Vector{MOI.VariableIndex},Vector{Vector{MOI.VariableIndex}}},
    }
    lagrangian_constraints::Vector{UMCT}
    constraint::MOI.ConstraintIndex{
        F,
        SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT},
    }
    monomials::MVT
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{
        SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT},
    },
    model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction,
    set::SOS.SOSPolynomialSet{<:SemialgebraicSets.BasicSemialgebraicSet},
) where {T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT}
    @assert MOI.output_dimension(f) == length(set.basis)
    # MOI does not modify the coefficients of the functions so we can modify `p`.
    # without altering `f`.
    # The monomials may be copied by MA however so we need to copy it.
    p = MB.algebra_element(MOI.Utilities.scalarize(f), copy(set.basis))
    λ_bases = B[]
    λ_variables =
        Union{Vector{MOI.VariableIndex},Vector{Vector{MOI.VariableIndex}}}[]
    λ_constraints = UMCT[]
    preprocessed =
        Certificate.preprocessed_domain(set.certificate, set.domain, p)
    for index in Certificate.preorder_indices(set.certificate, preprocessed)
        λ, λ_variable, λ_constraint, λ_basis = lagrangian_multiplier(
            model,
            set.certificate,
            index,
            preprocessed,
            T,
        )
        push!(λ_variables, λ_variable)
        push!(λ_constraints, λ_constraint)
        push!(λ_bases, λ_basis)
        # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
        # need to call `similar`. This is critical since `T` is
        # `Float64` when used with JuMP and the coefficient type is often `Int` if
        # `set.domain.V` is `FullSpace` or `FixedPolynomialSet`.
        g = Certificate.generator(set.certificate, index, preprocessed)
        # TODO replace with `MA.sub_mul` when it works.
        p = MA.operate!!(MA.add_mul, p, -one(T), λ, similar(g, T))
    end
    new_set = SOS.SOSPolynomialSet(
        set.domain.V,
        # For terms, `monomials` is `OneOrZeroElementVector`
        # so we convert it with `monomial_vector`
        # Later, we'll use `MP.MonomialBasis` which is going to do that anyway
        MB.SubBasis{MB.Monomial}(MP.monomial_vector(MP.monomials(p))),
        Certificate.ideal_certificate(set.certificate),
    )
    constraint = MOI.add_constraint(
        model,
        MOI.Utilities.vectorize(MP.coefficients(p)),
        new_set,
    )

    return SOSPolynomialInSemialgebraicSetBridge{
        T,
        F,
        DT,
        CT,
        B,
        UMCT,
        UMST,
        MCT,
        MT,
        MVT,
    }(
        λ_bases,
        λ_variables,
        λ_constraints,
        constraint,
        set.basis.monomials,
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
    ::Type{<:SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT}},
) where {T,F,DT,CT}
    return constrained_variable_types(SOS.matrix_cone_type(CT))
end
function MOI.Bridges.added_constraint_types(
    ::Type{
        SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT},
    },
) where {T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT}
    return [(F, SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT})]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialInSemialgebraicSetBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{
        <:SOS.SOSPolynomialSet{
            SemialgebraicSets.BasicSemialgebraicSet{S,PS,AT},
            MB.SubBasis{MB.Monomial,MT,MVT},
            CT,
        },
    },
) where {T,S,PS,AT,CT,MT,MVT}

    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOI.Utilities.promote_operation(-, T, F, MOI.VectorOfVariables)
    MCT = SOS.matrix_cone_type(CT)
    B = Certificate.multiplier_basis_type(CT, MT)
    UMCT = union_constraint_types(MCT)
    UMST = union_set_types(MCT)
    IC = Certificate.ideal_certificate(CT)
    return SOSPolynomialInSemialgebraicSetBridge{
        T,
        G,
        AT,
        IC,
        B,
        UMCT,
        UMST,
        MCT,
        MT,
        MVT,
    }
end

# Attributes, Bridge acting as an model
function MOI.get(
    bridge::SOSPolynomialInSemialgebraicSetBridge,
    ::MOI.NumberOfVariables,
)
    return mapreduce(_num_variables, +, bridge.lagrangian_variables, init = 0)
end
function MOI.get(
    bridge::SOSPolynomialInSemialgebraicSetBridge,
    ::MOI.ListOfVariableIndices,
)
    return Iterators.flatten([
        _list_variables(Qi) for Qi in bridge.lagrangian_variables
    ])
end
function MOI.get(
    bridge::SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,B,UMCT,UMST},
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
) where {T,F,DT,CT,B,UMCT,UMST,S<:UMST}
    return mapreduce(
        cQ ->
            _num_constraints(cQ, MOI.ConstraintIndex{MOI.VectorOfVariables,S}),
        +,
        bridge.lagrangian_constraints,
        init = 0,
    )
end
function MOI.get(
    b::SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,B,UMCT,UMST},
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {T,F,DT,CT,B,UMCT,UMST,S<:UMST}
    C = MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    return Iterators.flatten([
        _list_constraints(cQ, C) for cQ in bridge.lagrangian_constraints
    ])
end
function MOI.get(
    ::SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT},
    ::MOI.NumberOfConstraints{
        F,
        SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT},
    },
) where {T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT}
    return 1
end
function MOI.get(
    b::SOSPolynomialInSemialgebraicSetBridge{T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT},
    ::MOI.ListOfConstraintIndices{
        F,
        SOS.SOSPolynomialSet{DT,MB.SubBasis{MB.Monomial,MT,MVT},CT},
    },
) where {T,F,DT,CT,B,UMCT,UMST,MCT,MT,MVT}
    return [b.constraint]
end

# Indices
function MOI.delete(
    model::MOI.ModelLike,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    MOI.delete(model, bridge.constraint)
    for variables in bridge.lagrangian_variables
        _delete_variables(model, variables)
    end
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
    throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    dual = MOI.get(model, attr, bridge.constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.constraint)
    μ = MultivariateMoments.moment_vector(dual, set.basis)
    return [dot(mono, μ) for mono in bridge.monomials]
end
function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    return MOI.get(model, attr, bridge.constraint)
end

function MOI.get(
    model::MOI.ModelLike,
    attr::Union{
        SOS.CertificateBasis,
        SOS.GramMatrixAttribute,
        SOS.MomentMatrixAttribute,
    },
    bridge::SOSPolynomialInSemialgebraicSetBridge,
)
    return MOI.get(model, attr, bridge.constraint)
end
function _gram(
    f::Function,
    Qs::Vector{Vector{MOI.VariableIndex}},
    gram_bases,
    T::Type,
    MCT,
)
    return SOS.build_gram_matrix(gram_bases, MCT, T) do i
        return convert(Vector{T}, f(Qs[i]))
    end
end
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.LagrangianMultipliers,
    bridge::SOSPolynomialInSemialgebraicSetBridge{T,F,CT,B,UMCT,UMST,MCT,M},
) where {T,F,CT,B,UMCT,UMST,MCT,M}
    @assert eachindex(bridge.lagrangian_variables) ==
            eachindex(bridge.lagrangian_bases)
    return map(
        i -> _gram(
            Q -> MOI.get(model, MOI.VariablePrimal(attr.result_index), Q),
            bridge.lagrangian_variables[i],
            bridge.lagrangian_bases[i],
            T,
            MCT,
        ),
        eachindex(bridge.lagrangian_variables),
    )
end
