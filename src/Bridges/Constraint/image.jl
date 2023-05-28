struct ImageBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    UMCT<:Union{
        Vector{<:MOI.ConstraintIndex{MOI.VectorOfVariables}},
        MOI.ConstraintIndex{MOI.VectorOfVariables},
    },
    UMST,
    MCT,
    GB<:Union{
        Vector{<:Vector{<:MultivariateBases.AbstractPolynomialBasis}}, # Symmetry
        Vector{<:MultivariateBases.AbstractPolynomialBasis},           # Sparsity
        MultivariateBases.AbstractPolynomialBasis,                     # No reduction
    },
    ZB<:MultivariateBases.AbstractPolynomialBasis,
    CT<:SOS.Certificate.AbstractIdealCertificate,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
} <: MOI.Bridges.Constraint.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    constraint::MOI.ConstraintIndex{F}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{GeometricBridge{T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}},
    model::MOI.ModelLike,
    g::MOI.AbstractVectorFunction,
    s::SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace},
) where {T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}
    @assert MOI.output_dimension(g) == length(s.monomials)
    scalars = MOI.Utilities.eachscalar(g)
    gram_basis = SOS.Certificate.gram_basis(
        s.certificate,
        SOS.Certificate.with_variables(p, s.domain),
    )
    set = SOS.matrix_cone(MCT, length(gram_basis))
    f = MOI.Utilities.zero_with_output_dimension(F, MOI.dimension(set))
    k = 0
    found = Dict{eltype(gram_basis.monomials),Int}()
    used = falses(length(scalars))
    for j in eachindex(gram_basis.monomials)
        for i in 1:j
            k += 1
            mono = gram_basis.monomials[i] * gram_basis.monomials[j]
            if haskey(found, mono)
                var = MOI.add_variable(model)
                push!(variables, var)
                MOI.Utilities.operate_output_index!(-, T, found[mono], f, var)
                MOI.Utilities.operate_output_index!(+, T, k, f, var)
            else
                found[mono] = k
                t = MP.searchsortedfirst(s.monomials, mono, rev = true)
                if t in eachindex(s.monomials) && s.monomials[t] == mono
                    used[t] = true
                    MOI.Utilities.operate_output_index!(+, T, k, f, scalars[t])
                end
            end
        end
    end
    for i in eachindex(scalars)
        if !used[i]
            error("Infeasible")
        end
    end
    constraint = MOI.add_constraint(model, f, set)
    return GeometricBridge{T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}(
        variables,
        constraint,
        gram_basis,
        zero_constraint,
        s.domain,
        s.monomials,
        s.certificate,
    )
end

function MOI.supports_constraint(
    ::Type{GeometricBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,M,MV,C}},
) where {T}
    return Certificate.gram_basis_type <: MB.MonomialBasis
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:GeometricBridge{T,F,DT,UMCT,UMST,MCT}},
) where {T,F,DT,UMCT,UMST,MCT}
    return [MOI.Reals]
end

function MOI.Bridges.added_constraint_types(
    ::Type{GeometricBridge{T,F,G,CT}},
) where {T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}
    return [(F, CT)]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:GeometricBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}},
) where {T,MT,MVT,CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    F = MOI.Utilities.promote_operation(-, T, G, MOI.VectorOfVariables)
    MCT = SOS.matrix_cone_type(CT)
    UMST = union_set_types(MCT)
    GB = SOS.Certificate.gram_basis_type(CT)
    return GeometricBridge{T,G,DT,UMCT,UMST,MCT,GB}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::GeometricBridge, ::MOI.NumberOfVariables)
    return length(bridge.variables)
end
function MOI.get(bridge::GeometricBridge, ::MOI.ListOfVariableIndices)
    return bridge.variables
end
function MOI.get(
    bridge::GeometricBridge{T,F},
    ::MOI.NumberOfConstraints{F,S},
) where {T,F,S}
    return bridge.constraint isa MOI.ConstraintIndex{F,S} ? 1 : 0
end

function MOI.get(
    bridge::GeometricBridge{T,F},
    ::MOI.ListOfConstraintIndices{F,S},
) where {T,F,S}
    if bridge.constraint isa MOI.ConstraintIndex{F,S}
        return [bridge.constraint]
    else
        return MOI.ConstraintIndex{F,S}[]
    end
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::GeometricBridge)
    MOI.delete(model, bridge.constraint)
    if !isempty(bridge.variables)
        MOI.delete(model, bridge.variables)
    end
    return
end

# Attributes, Bridge acting as a constraint
function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::GeometricBridge,
)
    return bridge.set
end
function MOI.get(::MOI.ModelLike, ::MOI.ConstraintPrimal, ::GeometricBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::GeometricBridge{T},
) where {T}
    dual = MOI.get(model, attr, bridge.constraint)
    output = similar(dual, length(bridge.set.monomials))
    for i in eachindex(bridge.set.monomials)
        output[i] = dual[bridge.first[i]]
    end
    return output
end

function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::GeometricBridge,
)
    return MOI.get(model, MOI.ConstraintDual(attr.N), bridge.constraint)
end

function MOI.get(
    ::MOI.ModelLike,
    ::SOS.CertificateBasis,
    bridge::GeometricBridge,
)
    return bridge.gram_basis
end

function _gram(
    f::Function,
    Q::Vector{MOI.VariableIndex},
    gram_basis,
    T::Type,
    MCT,
)
    return SOS.build_gram_matrix(convert(Vector{T}, f(Q)), gram_basis, MCT, T)
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
    attr::SOS.GramMatrixAttribute,
    bridge::GeometricBridge{T,F,DT,UMCT,UMST,MCT},
) where {T,F,DT,UMCT,UMST,MCT}
    q = MOI.get(model, MOI.ConstraintPrimal(), model.constraint)
    return SOS.build_gram_matrix(q, model.gram_basis, MCT, T)
end
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.MomentMatrixAttribute,
    bridge::GeometricBridge,
)
    if bridge.cQ isa Vector{<:MOI.ConstraintIndex}
        return SOS.build_moment_matrix(bridge.gram_basis) do i
            return MOI.get(model, MOI.ConstraintDual(attr.N), bridge.cQ[i])
        end
    else
        return SOS.build_moment_matrix(
            MOI.get(model, MOI.ConstraintDual(attr.N), bridge.cQ),
            bridge.gram_basis,
        )
    end
end
