struct ScaledDiagonallyDominantBridge{T, F} <: MOIB.AbstractBridge
    side_dimension::Int
    equality::Vector{MOI.ConstraintIndex{F, MOI.EqualTo{T}}}
    psd2x2::Vector{MOI.ConstraintIndex{MOI.VectorAffineFunction{T},
                                       PositiveSemidefinite2x2ConeTriangle}}
end

function ScaledDiagonallyDominantBridge{T, F}(
    model::MOI.ModelLike, f::MOI.AbstractVectorFunction,
    s::ScaledDiagonallyDominantConeTriangle) where {T, F}
    # `p.Q` is SDD iff it is the sum of psd matrices Mij that are zero except for
    # entries ii, ij and jj [Lemma 9, AM17].
    @assert MOI.output_dimension(f) == MOI.dimension(s)
    @assert s.side_dimension > 1 # The bridges does not work with 1,
                                 # `matrix_add_constraint` should have used
                                 # Nonnegatives instead
    n = s.side_dimension
    fs = MOIU.eachscalar(f)
    # `g[r, c]` will contain the expression `f[r, c] - sum Mij[r, c]`
    # Cannot use `collect(fs)` as its `eltype` might be different to `F`, e.g.
    # if `f` is a `MOI.VectorOfVariables`.
    g = F[zero(F) for i in 1:MOI.dimension(s)]
    psd2x2 = Vector{MOI.ConstraintIndex{
        MOI.VectorAffineFunction{T}, PositiveSemidefinite2x2ConeTriangle}}(
            undef, length(fs) - n)
    diag_idx(i) = div(i * (i + 1), 2)
    k = 0
    k2x2 = 0
    for j in 1:n
        for i in 1:(j-1)
            k += 1
            MOIU.operate!(+, T, g[k], fs[k])
            k2x2 += 1
            Mii = MOI.SingleVariable(MOI.add_variable(model))
            MOIU.operate!(-, T, g[diag_idx(i)], Mii)
            Mij = MOI.SingleVariable(MOI.add_variable(model))
            MOIU.operate!(-, T, g[k], Mij)
            Mjj = MOI.SingleVariable(MOI.add_variable(model))
            MOIU.operate!(-, T, g[diag_idx(j)], Mjj)
            # PSD constraints on 2x2 matrices are SOC representable
            psd2x2[k2x2] = psd2x2_constraint(model, Mii, Mij, Mjj)
        end
        k += 1
        MOIU.operate!(+, T, g[k], fs[k])
    end
    equality = map(f -> MOIU.add_scalar_constraint(model, f, MOI.EqualTo(0.0)), g)
    return ScaledDiagonallyDominantBridge{T, F}(n, equality, psd2x2)
end

function MOI.supports_constraint(::Type{<:ScaledDiagonallyDominantBridge},
                                 ::Type{<:MOI.AbstractVectorFunction},
                                 ::Type{<:ScaledDiagonallyDominantConeTriangle})
    return true
end
function MOIB.added_constraint_types(::Type{ScaledDiagonallyDominantBridge{T, F}}) where {T, F}
    return [(F, MOI.EqualTo{T}),
            (MOI.VectorAffineFunction{T},
             PositiveSemidefinite2x2ConeTriangle)]
end
function MOIB.concrete_bridge_type(::Type{<:ScaledDiagonallyDominantBridge{T}},
                                   F::Type{<:MOI.AbstractVectorFunction},
                                   ::Type{ScaledDiagonallyDominantConeTriangle}) where T
    S = MOIU.scalar_type(F)
    G = MOIU.promote_operation(-, T, S, MOI.SingleVariable)
    return ScaledDiagonallyDominantBridge{T, G}
end

# Attributes, Bridge acting as an model
function MOI.get(bridge::ScaledDiagonallyDominantBridge{T, F},
                 ::MOI.NumberOfConstraints{F, MOI.EqualTo{T}}) where {T, F}
    return length(bridge.equality)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge{T},
                 ::MOI.NumberOfConstraints{MOI.VectorAffineFunction{T},
                                           PositiveSemidefinite2x2ConeTriangle}) where {T, F}
    return length(bridge.psd2x2)
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge{T, F},
                 ::MOI.ListOfConstraintIndices{F, MOI.EqualTo{T}}) where {T, F}
    return bridge.equality
end
function MOI.get(bridge::ScaledDiagonallyDominantBridge{T},
                 ::MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{T},
                                               MOI.RotatedSecondOrderCone}) where {T}
    return bridge.psd2x2
end

# Indices
function MOI.delete(model::MOI.ModelLike, bridge::ScaledDiagonallyDominantBridge)
    for ci in bridge.equality
        MOI.delete(model, ci)
    end
    for ci in bridge.psd2x2
        MOI.delete(model, ci)
    end
    # TODO delete variables
end

# TODO ConstraintPrimal and ConstraintDual
function MOI.get(model::MOI.ModelLike, attr::MOI.ConstraintDual,
                 bridge::ScaledDiagonallyDominantBridge)
    dual = copy(MOI.get(model, attr, bridge.equality))
    # Need to divide by 2 because of the custom scalar product for this cone
    k = 0
    for j in 1:bridge.side_dimension
        for i in 1:(j-1)
            k += 1
            dual[k] /= 2
        end
        k += 1
    end
    return dual
end
