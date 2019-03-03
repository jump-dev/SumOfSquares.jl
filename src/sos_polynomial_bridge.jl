struct SOSPolynomialBridge{T, F <: MOI.AbstractVectorFunction,
                           DT <: AbstractSemialgebraicSet,
                           VBS <: AbstractVariableBridge, MCT,
                           BT <: PolyJuMP.AbstractPolynomialBasis,
                           MT <: MP.AbstractMonomial,
                           MVT <: AbstractVector{MT}} <: MOIB.AbstractBridge
    variable_bridge::VBS
    certificate_monomials::MVT
    zero_constraint::MOI.ConstraintIndex{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}
end

function SOSPolynomialBridge{T, F, DT, VBS, MCT, BT, MT, MVT}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::SOSPolynomialSet{<:AbstractAlgebraicSet}) where {
        # Need to specify types to avoid ambiguity with the default constructor
        T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
        VBS <: AbstractVariableBridge, MCT,
        BT <: PolyJuMP.AbstractPolynomialBasis, MT <: MP.AbstractMonomial,
        MVT <: AbstractVector{MT}
    }
    @assert MOI.output_dimension(f) == length(s.monomials)
    p = MP.polynomial(collect(MOIU.eachscalar(f)), s.monomials)
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = convert(typeof(p), rem(p, ideal(s.domain)))
    X = monomials_half_newton_polytope(MP.monomials(r), s.newton_polytope)
    Q, variable_bridge = add_matrix_variable_bridge(model, MCT, length(X), T)
    g = build_gram_matrix(Q, X)
    q = r - g
    set = PolyJuMP.ZeroPolynomialSet(s.domain, s.basis, MP.monomials(q))
    coefs = MOIU.vectorize(MP.coefficients(q))
    zero_constraint = MOI.add_constraint(model, coefs, set)
    return SOSPolynomialBridge{T, F, DT, VBS, MCT, BT, MT, MVT}(
        variable_bridge, X, zero_constraint)
end

function MOI.supports_constraint(::Type{SOSPolynomialBridge{T}},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:SOSPolynomialSet{<:AbstractAlgebraicSet}}) where T
    return true
end
function MOIB.added_constraint_types(::Type{SOSPolynomialBridge{T, F, DT, VBS, MCT, BT, MT, MVT}}) where {T, F, DT, VBS, MCT, BT, MT, MVT}
    added = [(F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT})]
    return append_added_constraint_types(added, MCT, T)
end
function MOIB.concrete_bridge_type(::Type{<:SOSPolynomialBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{<:SOSPolynomialSet{DT, CT, <:PolyJuMP.MonomialBasis, MT, MVT}}) where {T, DT<:AbstractAlgebraicSet, CT, MT, MVT}
    # promotes VectorOfVariables into VectorAffineFunction, it should be enough
    # for most use cases
    G = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    MCT = matrix_cone_type(CT)
    VBS = union_vector_bridge_types(MCT, T)
    return SOSPolynomialBridge{T, G, DT, VBS, MCT, PolyJuMP.MonomialBasis, MT, MVT}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::SOSPolynomialBridge, attr::MOI.NumberOfVariables)
    return MOI.get(bridge.variable_bridge, attr)
end
function MOI.get(::SOSPolynomialBridge{T, F, DT, VBS, MCT, BT, MT, MVT},
                 ::MOI.NumberOfConstraints{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {
        # Need to specify types to avoid ambiguity with the method redirecting
        # to `variable_bridge`
        T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
        VBS <: AbstractVariableBridge, MCT <: MOI.AbstractVectorSet,
        BT <: PolyJuMP.AbstractPolynomialBasis, MT <: MP.AbstractMonomial,
        MVT <: AbstractVector{MT}
    }
    return 1
end
function MOI.get(b::SOSPolynomialBridge{T, F, DT, VBS, MCT, BT, MT, MVT},
                 ::MOI.ListOfConstraintIndices{F, PolyJuMP.ZeroPolynomialSet{DT, BT, MT, MVT}}) where {
        # Need to specify types to avoid ambiguity with the method redirecting
        # to `variable_bridge`
        T, F <: MOI.AbstractVectorFunction, DT <: AbstractSemialgebraicSet,
        VBS <: AbstractVariableBridge, MCT <: MOI.AbstractVectorSet,
        BT <: PolyJuMP.AbstractPolynomialBasis, MT <: MP.AbstractMonomial,
        MVT <: AbstractVector{MT}
    }
    return [b.zero_constraint]
end
function MOI.get(bridge::SOSPolynomialBridge,
                 attr::MOI.NumberOfConstraints)
    return MOI.get(bridge.variable_bridge, attr)
end
function MOI.get(bridge::SOSPolynomialBridge,
                 attr::MOI.ListOfConstraintIndices)
    return MOI.get(bridge.variable_bridge, attr)
end


# Indices
function MOI.delete(model::MOI.ModelLike, bridge::SOSPolynomialBridge)
    # First delete the constraints in which the Gram matrix appears
    MOI.delete(model, bridge.zero_constraint)
    # Now we delete the Gram matrix
    MOI.delete(model, bridge.variable_bridge)
end

struct ValueNotSupported <: Exception end
function Base.showerror(io::IO, ::ValueNotSupported)
    print(io, "`value` is no supported for Sum-of-Squares constraints, use",
          " `gram_matrix` instead.")
end

struct DualNotSupported <: Exception end
function Base.showerror(io::IO, ::DualNotSupported)
    print(io, "`dual` is no supported for Sum-of-Squares constraints in a",
          " domain, use `moment_matrix` instead.")
end


# Attributes, Bridge acting as a constraint
function MOI.get(::MOI.ModelLike,
                 ::MOI.ConstraintPrimal,
                 ::SOSPolynomialBridge)
    throw(ValueNotSupported())
end
function MOI.get(::MOI.ModelLike,
                 ::MOI.ConstraintDual,
                 ::SOSPolynomialBridge)
    throw(DualNotSupported())
end
function MOI.get(model::MOI.ModelLike,
                 attr::MOI.ConstraintDual,
                 bridge::SOSPolynomialBridge{T, <:MOI.AbstractVectorFunction, FullSpace}) where {T}
    return MOI.get(model, attr, bridge.zero_constraint)
end
function MOI.get(::MOI.ModelLike, ::CertificateMonomials,
                 bridge::SOSPolynomialBridge)
    return bridge.certificate_monomials
end
function MOI.get(model::MOI.ModelLike,
                 attr::GramMatrixAttribute,
                 bridge::SOSPolynomialBridge)
    return build_gram_matrix(MOI.get(model, attr, bridge.variable_bridge),
                             bridge.certificate_monomials)
end
function MOI.get(model::MOI.ModelLike,
                 attr::MomentMatrixAttribute,
                 bridge::SOSPolynomialBridge)
    return build_moment_matrix(MOI.get(model, attr, bridge.variable_bridge),
                               bridge.certificate_monomials)
end
