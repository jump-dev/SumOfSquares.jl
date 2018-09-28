export DSOSPoly, SDSOSPoly, SOSPoly

function JuMP.result_value(p::MatPolynomial{JuMP.VariableRef})
    MatPolynomial(map(JuMP.result_value, p.Q), p.x)
end

for poly in (:DSOSPoly, :SDSOSPoly, :SOSPoly)
    @eval begin
        struct $poly{PB<:PolyJuMP.AbstractPolynomialBasis} <: PolyJuMP.AbstractPoly
            polynomial_basis::PB
        end
        $poly(x::AbstractVector{<:MultivariatePolynomials.AbstractPolynomialLike}) = $poly(MonomialBasis(x))
    end
end

const PosPoly{PB} = Union{DSOSPoly{PB}, SDSOSPoly{PB}, SOSPoly{PB}}

JuMP.variable_type(m::JuMP.Model, p::PosPoly) = PolyJuMP.polytype(m, p, p.polynomial_basis)
PolyJuMP.polytype(m::JuMP.Model, ::PosPoly, basis::PolyJuMP.MonomialBasis{MT, MV}) where {MT<:AbstractMonomial, MV<:AbstractVector{MT}} = MatPolynomial{JuMP.VariableRef, MT, MV}

# Sum-of-Squares polynomial

_polytype(m::JuMP.Model, ::PosPoly, x::MVT) where {MT<:AbstractMonomial, MVT<:AbstractVector{MT}} = MatPolynomial{JuMP.VariableRef, MT, MVT}

"""
    constraint_matpoly!(model::JuMP.AbstractModel, p::MatPolynomial, ::Union{SOSPoly, SDSOSPoly, DSOSPoly})

Constraints matrix of `p` to be

* positive semidefinite if the third argument is of type `SOSPoly`,
* scaled diagonally dominant if the third argument is of type `SDSOSPoly`, or
* diagonally dominant if the third argument is of type `DSOSPoly`.

See Definition 4 of [AM17] for a precise definition of the last two items.

[AM17] Ahmadi, A. A. & Majumdar, A.
*DSOS and SDSOS Optimization: More Tractable Alternatives to Sum of Squares and Semidefinite Optimization*
ArXiv e-prints, 2017
"""
function constraint_matpoly! end

function constraint_matpoly!(m, p::MatPolynomial, ::SOSPoly)
    JuMP.add_constraint(m, JuMP.VectorConstraint(p.Q.Q, MOI.PositiveSemidefiniteConeTriangle(length(p.x))))
end
function constraint_matpoly!(model, p::MatPolynomial, ::SDSOSPoly)
    # `p.Q` is SDD iff it is the sum of psd matrices Mij that are zero except for
    # entries ii, ij and jj [Lemma 9, AM17].
    n = length(p.x)
    T = JuMP.GenericAffExpr{Float64, JuMP.VariableRef}
    # `Q[r, c]` will contain the expression `p.Q[r, c] - sum Mij[r, c]`
    Q = SymMatrix{T}(Vector{T}(p.Q.Q), n)
    for i in 1:n
        for j in i+1:n
            Mii = @variable(model)
            JuMP.add_to_expression!(Q[i, i], -1.0, Mii)
            Mij = @variable(model)
            JuMP.add_to_expression!(Q[i, j], -1.0, Mij)
            Mjj = @variable(model)
            JuMP.add_to_expression!(Q[j, j], -1.0, Mjj)
            # PSD constraints on 2x2 matrices are SOC representable
            @constraint(model, [Mii + Mjj, 2Mij, Mii - Mjj] in MOI.SecondOrderCone(3))
        end
    end
    @constraint(model, Q.Q .== 0)
end
function constraint_matpoly!(model, p::MatPolynomial, ::DSOSPoly)
    n = length(p.x)
    Q = Matrix{JuMP.VariableRef}(undef, n, n)
    for i in 1:n
        for j in i:n
            if i == j
                Q[i, j] = p[i, j]
            else
                Q[j, i] = Q[i, j] = JuMP.VariableRef(model)
                @constraint model Q[i, j] >= p[i, j]
                @constraint model Q[i, j] >= -p[i, j]
            end
        end
    end
    for i in 1:n
        @constraint model 2Q[i, i] >= sum(Q[i, :])
    end
    # If n > 1, this is implied by the constraint but it doesn't hurt to add the variable cone
    # Adding things on varCones makes JuMP think that it is SDP
    # push!(m.varCones, (:NonNeg, map(i -> p[i, i].col, 1:n)))
end
function _matpolynomial(m, x::AbstractVector{<:AbstractMonomial}, binary::Bool, integer::Bool)
    if isempty(x)
        zero(JuMP.variable_type(m, SOSPoly(x)))
    else
        function _newvar(i, j)
            v = JuMP.VariableRef(m)
            if length(x) == 1
                # 1x1 matrix is SDP iff its only entry is nonnegative
                # We handle this case here and do not create any SDP constraint
                setlowerbound(v, 0)
            end
            if binary
                setbinary(v)
            end
            if integer
                setinteger(v)
            end
            v
        end
        MatPolynomial{JuMP.VariableRef}(_newvar, x)
    end
end
function _createpoly(m::JuMP.Model, set::PosPoly, basis::PolyJuMP.MonomialBasis, binary::Bool, integer::Bool)
    p = _matpolynomial(m, basis.monomials, binary, integer)
    if length(basis.monomials) > 1
        constraint_matpoly!(m, p, set)
    end
    p
end
function PolyJuMP.createpoly(m::JuMP.Model, p::PosPoly, binary::Bool, integer::Bool)
    _createpoly(m, p, p.polynomial_basis, binary, integer)
end

# Defer other methods to defaults in PolyJuMP
createpoly(args...) = PolyJuMP.createpoly(args...)
