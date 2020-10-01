struct SOSPolynomialBridge{
    T, F <: MOI.AbstractVectorFunction,
    DT <: SemialgebraicSets.AbstractSemialgebraicSet,
    UMCT <: Union{Vector{<:MOI.ConstraintIndex{MOI.VectorOfVariables}},
                  MOI.ConstraintIndex{MOI.VectorOfVariables}},
    UMST, MCT,
    GB <: Union{Vector{<:MultivariateBases.AbstractPolynomialBasis},
                MultivariateBases.AbstractPolynomialBasis},
    ZB <: MultivariateBases.AbstractPolynomialBasis,
    CT <: SOS.Certificate.AbstractIdealCertificate,
    MT <: MP.AbstractMonomial,
    MVT <: AbstractVector{MT}} <: MOIB.Constraint.AbstractBridge

    Q::Union{Vector{Vector{MOI.VariableIndex}},
             Vector{MOI.VariableIndex}} # Vector{Vector{MOI.VariableIndex}} for sparse SOS
    cQ::UMCT
    gram_basis::GB
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, ZB, MT, MVT}}
    domain::DT
    monomials::MVT
    certificate::CT
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT}},
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}) where {T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT}

    @assert MOI.output_dimension(f) == length(s.monomials)
    # MOI does not modify the coefficients of the functions so we can modify `p`.
    # without altering `f`.
    # The monomials may be copied by MA however so we need to copy it.
    p = MP.polynomial(MOIU.scalarize(f), copy(s.monomials))
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
    # need to call `changecoefficienttype`. This is critical since `T` is
    # `Float64` when used with JuMP and the coefficient type is often `Int` if
    # `set.domain.V` is `FullSpace` or `FixedPolynomialsSet`.
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = SOS.Certificate.get(s.certificate, SOS.Certificate.ReducedPolynomial(), p, MP.changecoefficienttype(s.domain, T))
    gram_basis = SOS.Certificate.get(s.certificate, SOS.Certificate.GramBasis(), r)
    g, Q, cQ = SOS.add_gram_matrix(model, MCT, gram_basis, T)
    # MOI does not modify the coefficients of the functions so we can modify `r`.
    # without altering `f`.
    q = MA.operate!(-, r, g)
    set = PolyJuMP.ZeroPolynomialSet(s.domain, SOS.Certificate.zero_basis(s.certificate), MP.monomials(q))
    coefs = MOIU.vectorize(MP.coefficients(q))
    zero_constraint = MOI.add_constraint(model, coefs, set)
    return SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT}(
        Q, cQ, gram_basis, zero_constraint, s.domain, s.monomials, s.certificate)
end

function MOI.supports_constraint(::Type{SOSPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOS.SOSPolynomialSet{<:SemialgebraicSets.AbstractAlgebraicSet}}) where T
    return true
end
function MOIB.added_constrained_variable_types(::Type{<:SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT}}) where {T, F, DT, UMCT, UMST, MCT}
    return constrained_variable_types(MCT)
end
function MOIB.added_constraint_types(::Type{SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT}}) where {T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT}
    return [(F, PolyJuMP.ZeroPolynomialSet{DT, ZB, MT, MVT})]
end
function MOIB.Constraint.concrete_bridge_type(
    ::Type{<:SOSPolynomialBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:SOS.SOSPolynomialSet{DT, MT, MVT, CT}}) where {T, DT<:SemialgebraicSets.AbstractAlgebraicSet, MT, MVT, CT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    MCT = SOS.matrix_cone_type(CT)
    UMCT = union_constraint_types(MCT)
    UMST = union_set_types(MCT)
    GB = SOS.Certificate.get(CT, SOS.Certificate.GramBasisType())
    ZB = SOS.Certificate.zero_basis_type(CT)
    return SOSPolynomialBridge{T, G, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT}
end

# Attributes, Bridge acting as an model
_num_variables(Q::Vector{MOI.VariableIndex}) = length(Q)
_num_variables(Q::Vector{Vector{MOI.VariableIndex}}) = mapreduce(length, +, Q, init = 0)
function MOI.get(bridge::SOSPolynomialBridge, ::MOI.NumberOfVariables)
    _num_variables(bridge.Q)
end
_list_variables(Q::Vector{MOI.VariableIndex}) = Q
_list_variables(Q::Vector{Vector{MOI.VariableIndex}}) = Iterators.flatten(Q)
function MOI.get(bridge::SOSPolynomialBridge, ::MOI.ListOfVariableIndices)
    _list_variables(bridge.Q)
end
_num_constraints(cQ::Vector, ::Type{C}) where C = count(ci -> ci isa C, cQ)
_num_constraints(cQ::C, ::Type{C}) where C = 1
_num_constraints(cQ, ::Type) = 0
function MOI.get(bridge::SOSPolynomialBridge{T, F, DT, UMCT, UMST},
                 ::MOI.NumberOfConstraints{MOI.VectorOfVariables, S}) where {T, F, DT, UMCT, UMST, S<:UMST}
    return _num_constraints(bridge.cQ, MOI.ConstraintIndex{MOI.VectorOfVariables, S})
end
_list_constraints(cQ::Vector, ::Type{C}) where C = filter(ci -> ci isa C, cQ)
_list_constraints(cQ::C, ::Type{C}) where C = [cQ]
_list_constraints(cQ, C::Type) = C[]
function MOI.get(bridge::SOSPolynomialBridge{T, F, DT, UMCT, UMST},
                 ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, S}) where {T, F, DT, UMCT, UMST, S<:UMST}
    return _list_constraints(bridge.cQ, MOI.ConstraintIndex{MOI.VectorOfVariables, S})
end
function MOI.get(::SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, PolyJuMP.ZeroPolynomialSet{DT, ZB, MT, MVT}}) where {
        T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT
    }
    return 1
end
function MOI.get(b::SOSPolynomialBridge{T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, PolyJuMP.ZeroPolynomialSet{DT, ZB, MT, MVT}}) where {
        T, F, DT, UMCT, UMST, MCT, GB, ZB, CT, MT, MVT
    }
    return [b.zero_constraint]
end

# Indices
_delete_variables(model, Q::Vector{MOI.VariableIndex}) = MOI.delete(model, Q)
function _delete_variables(model, Qs::Vector{Vector{MOI.VariableIndex}})
    for Q in Qs
        MOI.delete(model, Q)
    end
end
function MOI.delete(model::MOI.ModelLike, bridge::SOSPolynomialBridge)
    # First delete the constraints in which the Gram matrix appears
    MOI.delete(model, bridge.zero_constraint)
    # Now we delete the Gram matrix
    _delete_variables(model, bridge.Q)
end

# Attributes, Bridge acting as a constraint
function MOI.get(model::MOI.ModelLike,
                 ::MOI.ConstraintSet,
                 bridge::SOSPolynomialBridge)
    return SOS.SOSPolynomialSet(bridge.domain, bridge.monomials, bridge.certificate)
end
function MOI.get(::MOI.ModelLike,
                 ::MOI.ConstraintPrimal,
                 ::SOSPolynomialBridge)
    throw(SOS.ValueNotSupported())
end

function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge{T}) where T
    dual = MOI.get(model, attr, bridge.zero_constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.zero_constraint)
    μ = MultivariateMoments.measure(dual, set.monomials)
    function reduced(mono)
        p = MP.polynomial(mono, T)
        domain = MP.changecoefficienttype(bridge.domain, T)
        return SOS.Certificate.get(
            bridge.certificate, SOS.Certificate.ReducedPolynomial(), p, domain)
    end
    return [dot(reduced(mono), μ) for mono in bridge.monomials]
end
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge{T, <:MOI.AbstractVectorFunction, SemialgebraicSets.FullSpace}) where {T}
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(model::MOI.ModelLike, attr::PolyJuMP.MomentsAttribute,
                 bridge::SOSPolynomialBridge)
    return MOI.get(model, attr, bridge.zero_constraint)
end

function MOI.get(::MOI.ModelLike, ::SOS.CertificateBasis,
                 bridge::SOSPolynomialBridge)
    return bridge.gram_basis
end
function _gram(f::Function, Q::Vector{MOI.VariableIndex}, gram_basis, T::Type)
    return SOS.build_gram_matrix(f(Q), gram_basis, T)
end
function _gram(f::Function, Qs::Vector{Vector{MOI.VariableIndex}}, gram_bases, T::Type)
    return SOS.SparseGramMatrix([_gram(f, Q, gram_basis, T) for (Q, gram_basis) in zip(Qs, gram_bases)])
end
function MOI.get(model::MOI.ModelLike,
                 attr::SOS.GramMatrixAttribute,
                 bridge::SOSPolynomialBridge{T}) where T
    return _gram(Q -> MOI.get(model, MOI.VariablePrimal(attr.N), Q),
                 bridge.Q, bridge.gram_basis, T::Type)
end
function MOI.get(model::MOI.ModelLike,
                 attr::SOS.MomentMatrixAttribute,
                 bridge::SOSPolynomialBridge)
    if bridge.cQ isa Vector{<:MOI.ConstraintIndex}
        return MultivariateMoments.SparseMomentMatrix([
            SOS.build_moment_matrix(MOI.get(model, MOI.ConstraintDual(attr.N), cQ), monos)
            for (cQ, monos) in zip(bridge.cQ, bridge.gram_basis)
        ])
    else
        return SOS.build_moment_matrix(MOI.get(model, MOI.ConstraintDual(attr.N),
                                               bridge.cQ),
                                       bridge.gram_basis)
    end
end
