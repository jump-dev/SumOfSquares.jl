struct GeometricBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    MT<:MP.AbstractMonomial,
    MVT<:AbstractVector{MT},
    CT<:SOS.Certificate.AbstractIdealCertificate,
} <: MOI.Bridges.Constraint.AbstractBridge
    variables::Vector{MOI.VariableIndex}
    constraint::MOI.ConstraintIndex{F}
    set::SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}
    gram_basis::MP.MonomialBasis{MT,MVT}
    first::Vector{Int}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{GeometricBridge{T,F,MT,MVT,CT}},
    model::MOI.ModelLike,
    g::MOI.AbstractVectorFunction,
    s::SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace},
) where {T,F,CT,MT,MVT}
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
    first = Union{Nothing,Tuple{Int,Int}}[nothing for _ in eachindex(scalars)]
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
                    first[t] = (i, j)
                    MOI.Utilities.operate_output_index!(+, T, k, f, scalars[t])
                end
            end
        end
    end
    for t in eachindex(scalars)
        if isnothing(first[t])
            error("Infeasible")
        end
    end
    constraint = MOI.add_constraint(model, f, set)
    return GeometricBridge{T,F,MT,MVT,CT}(
        variables,
        constraint,
        gram_basis,
        s,
        first,
    )
end

function MOI.supports_constraint(
    ::Type{GeometricBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,CT,MT,MVT}},
) where {T,CT,MT,MVT}
    return Certificate.gram_basis_type(CT) <: MB.MonomialBasis
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:GeometricBridge},
)
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(
    ::Type{GeometricBridge{T,F}},
) where {T,F}
    return Tuple{Type,Type}[(F, MOI.PositiveSemidefiniteConeTriangle)] # TODO
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:GeometricBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}},
) where {T,MT,MVT,CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    F = MOI.Utilities.promote_operation(-, T, G, MOI.VectorOfVariables)
    return GeometricBridge{T,F,CT,MT,MVT}
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
    attr::Union{MOI.ConstraintDual,PolyJuMP.MomentsAttribute},
    bridge::GeometricBridge{T},
) where {T}
    dual = MOI.get(model, MOI.ConstraintDual(attr.result_index), bridge.constraint)
    output = similar(dual, length(bridge.set.monomials))
    for i in eachindex(bridge.set.monomials)
        output[i] = dual[bridge.first[i]]
    end
    return output
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
    bridge::GeometricBridge{T},
) where {T,F,DT,UMCT,UMST,MCT}
    q = MOI.get(model, MOI.ConstraintPrimal(attr.result_index), bridge.constraint)
    return SOS.build_gram_matrix(q, model.gram_basis, MOI.PositiveSemidefiniteConeTriangle, T) # TODO
end
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.MomentMatrixAttribute,
    bridge::GeometricBridge,
)
    return SOS.build_moment_matrix(
        MOI.get(model, MOI.ConstraintDual(attr.result_index), bridge.constraint),
        bridge.gram_basis,
    )
end
