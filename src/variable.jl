export DSOSPoly, SDSOSPoly, SOSPoly

function JuMP.value(p::GramMatrix{JuMP.VariableRef})
    GramMatrix(map(JuMP.value, p.Q), p.x)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{PB<:PolyJuMP.AbstractPolynomialBasis} <: PolyJuMP.AbstractPoly
            polynomial_basis::PB
        end
        $poly(x::AbstractVector{<:MP.APL}) = $poly(MonomialBasis(x))
    end
end

matrix_cone_type(::SOSPoly) = MOI.PositiveSemidefiniteConeTriangle
matrix_cone_type(::DSOSPoly) = DiagonallyDominantConeTriangle
matrix_cone_type(::SDSOSPoly) = ScaledDiagonallyDominantConeTriangle

const PosPoly{PB} = Union{DSOSPoly{PB}, SDSOSPoly{PB}, SOSPoly{PB}}

JuMP.variable_type(m::JuMP.Model, p::PosPoly) = PolyJuMP.polytype(m, p, p.polynomial_basis)
gram_eltype(::Union{DSOSPoly, SOSPoly}) = JuMP.VariableRef
gram_eltype(::SDSOSPoly) = JuMP.AffExpr # affine because of the variable bridge
function PolyJuMP.polytype(m::JuMP.Model, cone::PosPoly,
                           basis::PolyJuMP.MonomialBasis{MT, MV}) where {MT<:MP.AbstractMonomial, MV<:AbstractVector{MT}}
    return GramMatrix{gram_eltype(cone), MT, MV}
end

# Sum-of-Squares polynomial

_polytype(m::JuMP.Model, ::PosPoly, x::MVT) where {MT<:MP.AbstractMonomial, MVT<:AbstractVector{MT}} = GramMatrix{JuMP.VariableRef, MT, MVT}

function moi_add_variable(model::MOI.ModelLike, set::MOI.AbstractVectorSet,
                          binary::Bool, integer::Bool)
    VB = variable_bridge_type(typeof(set), Float64)
    Q, variable_bridge = add_variable_bridge(VB, model, set)
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
                           v::PolyJuMP.Variable{<:PosPoly{<:PolyJuMP.MonomialBasis}},
                           name::String="")
    monos = v.p.polynomial_basis.monomials
    set = matrix_cone(matrix_cone_type(v.p), length(monos))
    Q = moi_add_variable(backend(model), set, v.binary, v.integer)
    F = JuMP.jump_function_type(model, eltype(Q))
    return build_gram_matrix(
        F[JuMP.jump_function(model, vi) for vi in Q], monos)
end
