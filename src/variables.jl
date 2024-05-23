export DSOSPoly, SDSOSPoly, SOSPoly

function PolyJuMP.bridges(::Type{<:PositiveSemidefinite2x2ConeTriangle})
    return [(Bridges.Variable.PositiveSemidefinite2x2Bridge, Float64)]
end
function PolyJuMP.bridges(::Type{<:ScaledDiagonallyDominantConeTriangle})
    return [(Bridges.Variable.ScaledDiagonallyDominantBridge, Float64)]
end
function PolyJuMP.bridges(::Type{<:CopositiveInnerCone})
    return [(Bridges.Variable.CopositiveInnerBridge, Float64)]
end

function JuMP.value(p::GramMatrix{<:JuMP.AbstractJuMPScalar})
    return GramMatrix(map(JuMP.value, p.Q), p.basis)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{B<:SA.ExplicitBasis} <: PolyJuMP.AbstractPoly
            basis::B
        end
        $poly(monos::AbstractVector{<:_APL}) = $poly(MB.SubBasis{MB.Monomial}(monos))
    end
end

matrix_cone_type(::SOSPoly) = MOI.PositiveSemidefiniteConeTriangle
matrix_cone_type(::DSOSPoly) = DiagonallyDominantConeTriangle
matrix_cone_type(::SDSOSPoly) = ScaledDiagonallyDominantConeTriangle

const PosPoly{PB} = Union{DSOSPoly{PB},SDSOSPoly{PB},SOSPoly{PB}}

function moi_add_variable(
    model::MOI.ModelLike,
    set::MOI.AbstractVectorSet,
    binary::Bool,
    integer::Bool,
)
    Q, _ = MOI.add_constrained_variables(model, set)
    if binary
        for q in Q
            MOI.add_constraint(model, q, MOI.ZeroOne())
        end
    end
    if integer
        for q in Q
            MOI.add_constraint(model, q, MOI.Integer())
        end
    end
    return Q
end

function JuMP.add_variable(
    model::JuMP.AbstractModel,
    v::PolyJuMP.Variable{<:PosPoly},
    ::String = "",
)
    MCT = matrix_cone_type(v.p)
    set = matrix_cone(MCT, length(v.p.basis))
    # FIXME There is no variable bridge mechanism yet:
    #       https://github.com/jump-dev/MathOptInterface.jl/issues/710
    #       so there is no equivalent to `BridgeableConstraint`.
    #       Yet, we may need constraint bridges here if it goes through
    #       the `generic_variable_bridge`.
    #       We don't need the `ScaledDiagonallyDominantBridge` since it does
    #       not use the `generic_variable_bridge`.
    JuMP.add_bridge(model, Bridges.Variable.PositiveSemidefinite2x2Bridge)
    JuMP.add_bridge(model, Bridges.Variable.ScaledDiagonallyDominantBridge)
    JuMP.add_bridge(model, Bridges.Variable.CopositiveInnerBridge)
    JuMP.add_bridge(model, Bridges.Constraint.EmptyBridge)
    JuMP.add_bridge(model, Bridges.Constraint.PositiveSemidefinite2x2Bridge)
    JuMP.add_bridge(model, Bridges.Constraint.DiagonallyDominantBridge)
    Q = moi_add_variable(backend(model), set, v.binary, v.integer)
    return build_gram_matrix(
        JuMP.VariableRef[JuMP.VariableRef(model, vi) for vi in Q],
        v.p.basis,
        MCT,
        Float64,
    )
end
