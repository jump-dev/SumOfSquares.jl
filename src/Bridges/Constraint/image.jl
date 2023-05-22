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
    cQ::UMCT
    gram_basis::GB
    variables::MOI.ConstraintIndex{
        F,
        PolyJuMP.ZeroPolynomialSet{DT,ZB,MT,MVT},
    }
    monomials::MVT
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialBridge{T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}},
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
    return SOSPolynomialBridge{T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}(
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
    ::Type{SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,M,MV,C}},
) where {T}
    return Certificate.gram_basis_type <: MB.MonomialBasis
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:SOSPolynomialBridge{T,F,DT,UMCT,UMST,MCT}},
) where {T,F,DT,UMCT,UMST,MCT}
    return [MOI.Reals]
end
function MOI.Bridges.added_constraint_types(
    ::Type{SOSPolynomialBridge{T,F,G,CT}},
) where {T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}
    return [(F, CT)]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialBridge{T}},
    G::Type{<:MOI.AbstractVectorFunction},
    ::Type{SOS.SOSPolynomialSet{SemialgebraicSets.FullSpace,MT,MVT,CT}},
) where {T,MT,MVT,CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    F = MOI.Utilities.promote_operation(-, T, G, MOI.VectorOfVariables)
    MCT = SOS.matrix_cone_type(CT)
    UMST = union_set_types(MCT)
    GB = SOS.Certificate.gram_basis_type(CT)
    return SOSPolynomialBridge{T,G,DT,UMCT,UMST,MCT,GB}
end

# Attributes, Bridge acting as an model
_num_variables(Q::Vector{MOI.VariableIndex}) = length(Q)
function _num_variables(Q::Vector{Vector{MOI.VariableIndex}})
    return mapreduce(length, +, Q, init = 0)
end
function MOI.get(bridge::SOSPolynomialBridge, ::MOI.NumberOfVariables)
    return _num_variables(bridge.Q)
end
_list_variables(Q::Vector{MOI.VariableIndex}) = Q
_list_variables(Q::Vector{Vector{MOI.VariableIndex}}) = Iterators.flatten(Q)
function MOI.get(bridge::SOSPolynomialBridge, ::MOI.ListOfVariableIndices)
    return _list_variables(bridge.Q)
end
_num_constraints(cQ::Vector, ::Type{C}) where {C} = count(ci -> ci isa C, cQ)
_num_constraints(cQ::C, ::Type{C}) where {C} = 1
_num_constraints(cQ, ::Type) = 0
function MOI.get(
    bridge::SOSPolynomialBridge{T,F,DT,UMCT,UMST},
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
) where {T,F,DT,UMCT,UMST,S<:UMST}
    return _num_constraints(
        bridge.cQ,
        MOI.ConstraintIndex{MOI.VectorOfVariables,S},
    )
end
_list_constraints(cQ::Vector, ::Type{C}) where {C} = filter(ci -> ci isa C, cQ)
_list_constraints(cQ::C, ::Type{C}) where {C} = [cQ]
_list_constraints(cQ, C::Type) = C[]
function MOI.get(
    bridge::SOSPolynomialBridge{T,F,DT,UMCT,UMST},
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {T,F,DT,UMCT,UMST,S<:UMST}
    return _list_constraints(
        bridge.cQ,
        MOI.ConstraintIndex{MOI.VectorOfVariables,S},
    )
end
function MOI.get(
    ::SOSPolynomialBridge{T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT},
    ::MOI.NumberOfConstraints{F,PolyJuMP.ZeroPolynomialSet{DT,ZB,MT,MVT}},
) where {T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}
    return 1
end
function MOI.get(
    b::SOSPolynomialBridge{T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT},
    ::MOI.ListOfConstraintIndices{F,PolyJuMP.ZeroPolynomialSet{DT,ZB,MT,MVT}},
) where {T,F,DT,UMCT,UMST,MCT,GB,ZB,CT,MT,MVT}
    return [b.zero_constraint]
end

# Indices
function _delete_variables(model, Q::Vector{MOI.VariableIndex})
    if !isempty(Q)
        # FIXME Since there is not variables in the list, we cannot
        # identify the `EmptyBridge` to delete
        MOI.delete(model, Q)
    end
end
function _delete_variables(model, Qs::Vector{Vector{MOI.VariableIndex}})
    for Q in Qs
        _delete_variables(model, Q)
    end
end
function MOI.delete(model::MOI.ModelLike, bridge::SOSPolynomialBridge)
    # First delete the constraints in which the Gram matrix appears
    MOI.delete(model, bridge.zero_constraint)
    # Now we delete the Gram matrix
    return _delete_variables(model, bridge.Q)
end

# Attributes, Bridge acting as a constraint
function MOI.get(
    model::MOI.ModelLike,
    ::MOI.ConstraintSet,
    bridge::SOSPolynomialBridge,
)
    return SOS.SOSPolynomialSet(
        bridge.domain,
        bridge.monomials,
        bridge.certificate,
    )
end
function MOI.get(::MOI.ModelLike, ::MOI.ConstraintPrimal, ::SOSPolynomialBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::SOSPolynomialBridge{T},
) where {T}
    dual = MOI.get(model, attr, bridge.zero_constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.zero_constraint)
    μ = MultivariateMoments.measure(dual, set.monomials)
    function reduced(mono)
        p = MP.polynomial(mono, T)
        domain = MP.changecoefficienttype(bridge.domain, T)
        return SOS.Certificate.reduced_polynomial(bridge.certificate, p, domain)
    end
    return [dot(reduced(mono), μ) for mono in bridge.monomials]
end
function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::SOSPolynomialBridge{
        T,
        <:MOI.AbstractVectorFunction,
        SemialgebraicSets.FullSpace,
    },
) where {T}
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::SOSPolynomialBridge,
)
    return MOI.get(model, attr, bridge.zero_constraint)
end

function MOI.get(
    ::MOI.ModelLike,
    ::SOS.CertificateBasis,
    bridge::SOSPolynomialBridge,
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
    bridge::SOSPolynomialBridge{T,F,DT,UMCT,UMST,MCT},
) where {T,F,DT,UMCT,UMST,MCT}
    return _gram(
        Q -> MOI.get(model, MOI.VariablePrimal(attr.N), Q),
        bridge.Q,
        bridge.gram_basis,
        T::Type,
        MCT,
    )
end
function MOI.get(
    model::MOI.ModelLike,
    attr::SOS.MomentMatrixAttribute,
    bridge::SOSPolynomialBridge,
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
