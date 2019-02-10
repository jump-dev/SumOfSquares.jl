export DSOSPoly, SDSOSPoly, SOSPoly

function JuMP.value(p::MatPolynomial{JuMP.VariableRef})
    MatPolynomial(map(JuMP.value, p.Q), p.x)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{PB<:PolyJuMP.AbstractPolynomialBasis} <: PolyJuMP.AbstractPoly
            polynomial_basis::PB
        end
        $poly(x::AbstractVector{<:MultivariatePolynomials.AbstractPolynomialLike}) = $poly(MonomialBasis(x))
    end
end

matrix_cone(::SOSPoly) = MOI.PositiveSemidefiniteConeTriangle
matrix_cone(::DSOSPoly) = DiagonallyDominantConeTriangle
matrix_cone(::SDSOSPoly) = ScaledDiagonallyDominantConeTriangle

const PosPoly{PB} = Union{DSOSPoly{PB}, SDSOSPoly{PB}, SOSPoly{PB}}

JuMP.variable_type(m::JuMP.Model, p::PosPoly) = PolyJuMP.polytype(m, p, p.polynomial_basis)
PolyJuMP.polytype(m::JuMP.Model, ::PosPoly, basis::PolyJuMP.MonomialBasis{MT, MV}) where {MT<:AbstractMonomial, MV<:AbstractVector{MT}} = MatPolynomial{JuMP.VariableRef, MT, MV}

# Sum-of-Squares polynomial

_polytype(m::JuMP.Model, ::PosPoly, x::MVT) where {MT<:AbstractMonomial, MVT<:AbstractVector{MT}} = MatPolynomial{JuMP.VariableRef, MT, MVT}

"""
    matrix_polynomial(model::Union{JuMP.AbstractModel, MOI.ModelLike},
                      cone::PosPoly, basis::PolyJuMP.AbstractPolynomialBasis)

Returns a polynomial of the form ``x^\\top Q x`` where `x` is a vector of the
elements of the polynomial basis and `Q` is a symmetric matrix of `model`
variables constrained to belong to a subset of the cone of symmetric positive
semidefinite matrix determined by `PosPoly`.
"""
function matrix_polynomial(var_type::Type, new_var, monos, matrix_cone)
    p = MatPolynomial{var_type}(new_var, monos)

    return p
end

function JuMP.add_variable(model::JuMP.AbstractModel,
                           v::PolyJuMP.Variable{<:PosPoly{<:PolyJuMP.MonomialBasis}},
                           name::String="")
    monos = v.p.polynomial_basis.monomials
    function new_var(i, j)
        vref = JuMP.VariableRef(model)
        if v.binary
            JuMP.set_binary(vref)
        end
        if v.integer
            JuMP.set_integer(vref)
        end
        return vref
    end
    p = MatPolynomial{JuMP.VariableRef}(new_var, monos)
    matrix_add_constraint(model, p, matrix_cone(v.p))
    return p
end
