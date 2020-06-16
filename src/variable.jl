export DSOSPoly, SDSOSPoly, SOSPoly

function PolyJuMP.bridges(::Type{<:PositiveSemidefinite2x2ConeTriangle})
    return [Bridges.Variable.PositiveSemidefinite2x2Bridge]
end
function PolyJuMP.bridges(::Type{<:ScaledDiagonallyDominantConeTriangle})
    return [Bridges.Variable.ScaledDiagonallyDominantBridge]
end
function PolyJuMP.bridges(::Type{<:CopositiveInnerCone})
    return [Bridges.Variable.CopositiveInnerBridge]
end

function JuMP.value(p::GramMatrix{<:JuMP.AbstractJuMPScalar})
    GramMatrix(map(JuMP.value, p.Q), p.basis)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{PB<:AbstractPolynomialBasis} <: PolyJuMP.AbstractPoly
            polynomial_basis::PB
        end
        $poly(x::AbstractVector{<:MP.APL}) = $poly(MonomialBasis(x))
    end
end

matrix_cone_type(::SOSPoly) = MOI.PositiveSemidefiniteConeTriangle
matrix_cone_type(::DSOSPoly) = DiagonallyDominantConeTriangle
matrix_cone_type(::SDSOSPoly) = ScaledDiagonallyDominantConeTriangle

const PosPoly{PB} = Union{DSOSPoly{PB}, SDSOSPoly{PB}, SOSPoly{PB}}

function moi_add_variable(model::MOI.ModelLike, set::MOI.AbstractVectorSet,
                          binary::Bool, integer::Bool)
    Q, con_Q = MOI.add_constrained_variables(model, set)
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

function JuMP.add_variable(model::JuMP.AbstractModel,
                           v::PolyJuMP.Variable{<:PosPoly},
                           name::String="")
    set = matrix_cone(matrix_cone_type(v.p), length(v.p.polynomial_basis))
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
        v.p.polynomial_basis
    )
end
